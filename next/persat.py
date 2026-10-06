"""
persat.py -- SAT search for PERIODIC tilings admitted by the depth-k hierarchical contact rule L_k.
Torus Z^N / (a1 Z x ... x aN Z) (optionally sheared: period vectors given as a matrix).
A periodic admissible tiling that is not fully de-substitutable is a genuine counterexample to
enforcement by the rule (independent of the pairwise certificate).
usage: python3 persat.py NAME k a b c [class]
"""
import sys, time, itertools, pickle
sys.path.insert(0, '/home/user/research'); sys.path.insert(0, '/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from chairN import apply, compose, inverse, vadd, vsub, vscale, identity, make_group
from polytile import normalize_cells, dissections, dissection_classes
from rigidity import Hier
from cert3 import log
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType

def torus_search(H, L, k, periods, want_models=1, timeout=None, forbid_hier=False):
    N = H.N; G = make_group(N)
    marks = sorted({mk for (_, _, mk, _) in L} | {mk for (_, _, _, mk) in L}) if k else [()]
    # restrict orientations to those occurring in L (as g_rel composed from identity)... use all of G
    gens = [g for g in G]
    def red(c): return tuple(c[i] % periods[i] for i in range(N))
    cells_list = sorted(H.cells)
    # placements: (g, t, mk) with t in torus; cubes in CELL units (not doubled)
    ones = tuple([1] * N)
    def tile_cells(g, t):   # tau-convention (as in Hier/PolyTile): cell = g(c) + tau + (g(1)-1)/2
        s_ = apply(g, ones); off = tuple((s_[i] - 1) // 2 for i in range(N))
        return [vadd(vadd(apply(g, c), t), off) for c in cells_list]
    var = {}; plist = []
    for g in gens:
        for t in itertools.product(*[range(p) for p in periods]):
            cs = [red(c) for c in tile_cells(g, t)]
            if len(set(cs)) < len(cs): return None   # period too small for the tile
            for mk in marks:
                var[(g, t, mk)] = len(plist) + 1; plist.append((g, t, mk, cs))
    cover = {}
    for v, (g, t, mk, cs) in enumerate(plist, start=1):
        for c in cs: cover.setdefault(c, []).append(v)
    # allowed neighbours per (placement, panel): from L (g_rel, t_rel, mk, mk2) in cell units
    byA = {}
    for (gr, tr, mk, mk2) in L: byA.setdefault(mk, []).append((gr, tr, mk2))
    dirs = []
    for i in range(N):
        for s in (1, -1):
            d = [0] * N; d[i] = s; dirs.append(tuple(d))
    panels = [(c, d) for c in cells_list for d in dirs if vadd(c, d) not in H.cells]
    # for each mark, for each panel, list of relative neighbours (gr,tr,mk2) covering the panel cell
    cov = {}
    for mk in marks:
        lst = []
        for (c, d) in panels:
            x = vadd(c, d)
            lst.append([(gr, tr, mk2) for (gr, tr, mk2) in byA.get(mk, []) if x in set(tile_cells(gr, tr))])
        cov[mk] = lst
    S = Cadical153(); top = len(plist)
    # exact cover
    for c, vs in cover.items():
        S.add_clause(vs)
        enc = CardEnc.atmost(lits=vs, bound=1, top_id=top, encoding=EncType.seqcounter)
        top = max(top, enc.nv)
        for cl in enc.clauses: S.add_clause(cl)
    # neighbour constraints
    ncl = 0
    for v, (g, t, mk, cs) in enumerate(plist, start=1):
        for pi, (c, d) in enumerate(panels):
            opts = []
            for (gr, tr, mk2) in cov[mk][pi]:
                g2 = compose(g, gr); t2 = red(vadd(t, apply(g, tr)))
                w = var.get((g2, t2, mk2))
                if w: opts.append(w)
            S.add_clause([-v] + opts); ncl += 1
    models = []
    while len(models) < want_models and S.solve():
        m = S.get_model(); sol = [plist[v - 1][:3] for v in range(1, len(plist) + 1) if m[v - 1] > 0]
        models.append(sol)
        S.add_clause([-var[s] for s in sol])
    return models

def desubstitute(H, tiles, k, periods, maxlevel=8):
    """de-substitution on the torus with tau-convention poses.  k>=1: level-2 supertiles are read from
    the marks (each tile's parent must have all 2^N children present with the right indices); then
    level m -> m+1 by exact grouping into D-parents.  Returns (level reached, ok)."""
    N = H.N
    def red(c): return tuple(c[i] % periods[i] for i in range(N))
    def kids(S, m):            # children (level m-1) of level-m supertile S=(G,Tau)
        G, Tau = S
        return [(compose(G, h), red(vadd(Tau, apply(G, tau)))) for (h, tau) in H.kids_maps(H.D, m)]
    def parent(child, i, m):    # level-m supertile having `child` as child i
        h, tau = H.kids_maps(H.D, m)[i]
        G = compose(child[0], inverse(h)); return (G, red(vsub(child[1], apply(G, tau))))
    tset = {(g, red(t), mk) for (g, t, mk) in tiles}
    if k >= 1:
        by_pos = {}
        for (a, b, m2) in tset: by_pos.setdefault((a, b), set()).add(m2[0])
        lvl = set()
        for (g, t, mk) in tset:
            S = parent((g, t), mk[0], 2)
            if not all(i in by_pos.get(kd, ()) for i, kd in enumerate(kids(S, 2))): return 2, False
            lvl.add(S)
        level = 2
    else:
        lvl = {(g, red(t)) for (g, t, mk) in tset}; level = 1
    while level < maxlevel and len(lvl) > 1:
        rem = set(lvl); parents = set()
        while rem:
            s_ = next(iter(rem)); found = None
            for i in range(len(H.D)):
                S = parent(s_, i, level + 1); K = kids(S, level + 1)
                if all(kd in rem for kd in K): found = (S, K); break
            if not found: return level + 1, False
            rem -= set(found[1]); parents.add(found[0])
        lvl = parents; level += 1
    return level, True

if __name__ == '__main__':
    from census3 import SHAPES
    name = sys.argv[1]; k = int(sys.argv[2]); periods = tuple(int(x) for x in sys.argv[3:6])
    ci = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    cells = normalize_cells(frozenset(SHAPES[name])); G = make_group(3)
    classes, _, _ = dissection_classes(cells, G); D = classes[ci]
    H = Hier(cells, D); P = H.adjacency_language(); L = H.language(k)
    log(f"### {name} class {ci}: |P|={len(P)} |L_{k}|={len(L)}; torus {periods}")
    t0 = time.time()
    models = torus_search(H, L, k, periods, want_models=3)
    if models is None: log("period too small"); sys.exit()
    log(f"  {len(models)} periodic tilings found ({time.time()-t0:.0f}s)")
    for mi, sol in enumerate(models):
        lv, ok = desubstitute(H, sol, k, periods)
        log(f"  model {mi}: {len(sol)} tiles; de-substitution {'ok' if ok else 'FAILS'} at level {lv}")
    pickle.dump(dict(name=name, k=k, periods=periods, ci=ci, models=models), open(f'persat_{name}_k{k}_{"x".join(map(str,periods))}_c{ci}.pkl', 'wb'))

# ---------------------------------------------------------------- inflation lemma (verification)
def inflate(H, tiles, k, periods):
    """D-expand every depth-k marked tile of a periodic L_k-tiling: new marks (i, old mark)[:k+1];
    returns depth-(k+1) tiles with doubled periods."""
    out = []
    for (g, t, mk) in tiles:
        for i, (h, tau) in enumerate(H.kids_maps(H.D, 2)):
            out.append((compose(g, h), vadd(vscale(2, t), apply(g, tau)), ((i,) + tuple(mk))[:k + 1]))
    return out, tuple(2 * p for p in periods)

def check_admissible(H, L, tiles, periods):
    """exact cover of the torus + every face contact (lifted) in L."""
    N = H.N
    def red(c): return tuple(c[i] % periods[i] for i in range(N))
    ones = tuple([1] * N); cells_list = sorted(H.cells)
    def tile_cells(g, t):
        s_ = apply(g, ones); off = tuple((s_[i] - 1) // 2 for i in range(N))
        return [vadd(vadd(apply(g, c), t), off) for c in cells_list]
    occ = {}
    for (g, t, mk) in tiles:
        for c in tile_cells(g, t):
            if red(c) in occ: return False, 'overlap'
            occ[red(c)] = (g, t, mk, c)
    vol = 1
    for p in periods: vol *= p
    if len(occ) != vol: return False, 'gap'
    dirs = []
    for i in range(N):
        for s in (1, -1):
            d = [0] * N; d[i] = s; dirs.append(tuple(d))
    bad = 0
    for (g, t, mk) in tiles:
        gi = inverse(g)
        for c in tile_cells(g, t):
            for d in dirs:
                x = vadd(c, d); g2, t2, mk2, c2 = occ[red(x)]
                if (g2, red(t2), mk2) == (g, red(t), mk): continue
                # lift t2 so that the neighbour covers x exactly: shift by the period multiple x - c2
                sh = vsub(x, c2)
                if any(sh[i] % periods[i] for i in range(N)): pass
                t2l = vadd(t2, sh)
                if (compose(gi, g2), apply(gi, vsub(t2l, t)), mk, mk2) not in L: bad += 1
    return bad == 0, bad

def desub_Z(H, tiles, k, periods, maxlevel=6):
    """de-substitution in Z^N of the periodic tiling (lifted over enough periods).  Follows the
    supertile containing the origin cell: level 1 -> 2 via marks (k>=1) or grouping, then level m -> m+1
    by searching a D-parent all of whose children are present.  Returns (first level without parent, or
    maxlevel+1 if all found, ok flag)."""
    N = H.N
    def kids(S, m):
        G, Tau = S
        return [(compose(G, h), vadd(Tau, apply(G, tau))) for (h, tau) in H.kids_maps(H.D, m)]
    def parent(child, i, m):
        h, tau = H.kids_maps(H.D, m)[i]
        G = compose(child[0], inverse(h)); return (G, vsub(child[1], apply(G, tau)))
    # lift: enough copies so that a level-maxlevel supertile around the origin is inside
    bb = max(max(c) for c in H.cells) + 1
    R = (2 ** (maxlevel - 1) * bb * 2) // min(periods) + 2
    lifted = set()
    for (g, t, mk) in tiles:
        for sh in itertools.product(range(-R, R + 1), repeat=N):
            lifted.add((g, vadd(t, tuple(sh[i] * periods[i] for i in range(N))), mk))
    pos = {}
    for (g, t, mk) in lifted: pos.setdefault((g, t), set()).add(mk)
    ones = tuple([1] * N); cells_list = sorted(H.cells)
    def tile_cells(g, t):
        s_ = apply(g, ones); off = tuple((s_[i] - 1) // 2 for i in range(N))
        return [vadd(vadd(apply(g, c), t), off) for c in cells_list]
    # tile at the origin
    origin = tuple([0] * N); cur = None
    for (g, t, mk) in lifted:
        if origin in tile_cells(g, t): cur = (g, t, mk); break
    level = 1
    present1 = {(g, t) for (g, t, mk) in lifted}
    if k >= 1:
        S = parent(cur[:2], cur[2][0], 2)
        if not all(i in pos.get(kd, ()) for i, kd in enumerate(kids(S, 2))): return 2, False
        # all level-2 supertiles (from marks)
        lvl = set(parent((g, t), mk[0], 2) for (g, t, mk) in lifted); cur = S; level = 2
    else:
        lvl = present1; cur = cur[:2]
    while level < maxlevel:
        found = None
        for i in range(len(H.D)):
            S = parent(cur, i, level + 1)
            if all(kd in lvl for kd in kids(S, level + 1)): found = S; break
        if found is None: return level + 1, False
        # build the full next level set (only supertiles whose children are all present)
        nxt = set(); seen = set()
        for s_ in lvl:
            for i in range(len(H.D)):
                S = parent(s_, i, level + 1)
                if S in seen: continue
                seen.add(S)
                if all(kd in lvl for kd in kids(S, level + 1)): nxt.add(S)
        lvl = nxt; cur = found; level += 1
    return level, True
