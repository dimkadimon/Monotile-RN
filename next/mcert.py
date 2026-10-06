"""
mcert.py -- certificate for TWO-LEVEL markings: a tile is marked by its frame AND its index inside its
parent supertile.  Marked pose: (g, t, i).  Rule at level 1: set of (g_rel, t_rel, i_A, i_B).

Multi-level scheme (no union/closure at the tile level):
   R_1 = hierarchical contact language of marked tiles
   R_{j+1} = relative poses of two level-j supertiles (shape 2^j T, unmarked) that are pairwise admissible,
             i.e. all contacts between their children lie in R_j.
   stop when R_{j+1} == 2 R_j (for j >= 2; both unmarked)  -- the structure is self-similar from level j on.
   Then the enclosure analysis (unique parent + presence) is run at every level 1..j.  If all pass, every
   tiling by the marked tile is uniquely hierarchical => strongly aperiodic (lattice-registered setting).
Level j objects are represented as PolyTile(2^(j-1) * cells, scaled dissection) with 'index' 0 when unmarked.
"""
import sys, time, itertools, pickle
sys.path.insert(0, '/home/user/research'); sys.path.insert(0, '/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from chairN import apply, compose, inverse, vadd, vsub, vscale, identity, det
from polytile import PolyTile, normalize_cells, scaled_cells
from cert3 import log

def level_tile(cells, dissection, j):
    """PolyTile for level-j supertiles (shape 2^(j-1) T) with the dissection scaled accordingly."""
    f = 2 ** (j - 1)
    cs = scaled_cells(cells, f) if f > 1 else cells
    d = []
    for (g, t) in dissection:
        sg = apply(g, tuple([1] * len(t)))
        d.append((g, tuple(f * t[i] + (f - 1) * (1 - sg[i]) // 2 for i in range(len(t)))))
    return PolyTile(cs, d, name=f'L{j}')

class MRules:
    """shape + rule R (set of (g, t, iA, iB)); marked poses (g, t, i); marked == False -> all indices 0."""
    def __init__(self, ch, R, marked):
        self.ch = ch; self.N = ch.N; self.marked = marked
        self.R = frozenset(R)
        self.k = len(ch.children)
        self.idx = range(self.k) if marked else [0]
        self.dirs = ch.dirs
        self._cubes = {}; self._halo = {}
        self.local = tuple(ch.local_cubes); self.local_set = frozenset(ch.local_cubes)
        self.panels = ch.panels; self.pidx = ch.panel_index
        # relative candidates per index of A: rp = (g, t, iB)
        self.rel_by_i = {i: sorted({(g, t, ib) for (g, t, ia, ib) in self.R if ia == i}) for i in self.idx}
        self.cover = {}
        for i in self.idx:
            cov = [[] for _ in self.panels]
            for rp in self.rel_by_i[i]:
                cb = self.cubes(rp)
                for kk, (c, d) in enumerate(self.panels):
                    if vadd(c, d) in cb: cov[kk].append(rp)
            self.cover[i] = cov
        self.origin = (ch.id, tuple([0] * self.N), 0)
        self.sup_for_child = []
        for (h, tau) in ch.children:
            G = inverse(h); self.sup_for_child.append((G, tuple(-x for x in apply(G, tau))))

    def cubes(self, pose):
        key = (pose[0], pose[1])
        r = self._cubes.get(key)
        if r is None:
            g, t = key; tt = vscale(2, t)
            r = frozenset(vadd(apply(g, c), tt) for c in self.local); self._cubes[key] = r
        return r
    def halo(self, pose):
        key = (pose[0], pose[1]); r = self._halo.get(key)
        if r is None:
            cb = self.cubes(pose); r = frozenset(vadd(c, d) for c in cb for d in self.dirs) - cb; self._halo[key] = r
        return r
    def rel(self, A, B):
        gi = inverse(A[0]); return (compose(gi, B[0]), apply(gi, vsub(B[1], A[1])))
    def pair_ok(self, A, B):
        ca = self.cubes(A); cb = self.cubes(B)
        if not ca.isdisjoint(cb): return False
        if self.halo(A).isdisjoint(cb): return True
        g, t = self.rel(A, B); return (g, t, A[2], B[2]) in self.R
    def absolute(self, A, rp):
        gA, tA, _ = A; g, t, ib = rp
        return (compose(gA, g), vadd(tA, apply(gA, t)), ib)
    def local_panel_index(self, A, c, d):
        gi = inverse(A[0]); return self.pidx[(apply(gi, vsub(c, vscale(2, A[1]))), apply(gi, d))]
    def candidates_for(self, A, c, d):
        return [self.absolute(A, rp) for rp in self.cover[A[2]][self.local_panel_index(A, c, d)]]
    def kids(self, S):
        G, T = S
        return [(compose(G, h), vadd(T, apply(G, tau)), (j if self.marked else 0)) for j, (h, tau) in enumerate(self.ch.children)]
    def supertiles_containing(self, A):
        """parents (c, S, kids) in which A is child c.  With markings, c must equal A's index."""
        out = []; seen = set()
        slots = [A[2]] if self.marked else range(self.k)
        for c in slots:
            G0, T0 = self.sup_for_child[c]
            S = (compose(A[0], G0), vadd(A[1], apply(A[0], T0)))
            if S in seen: continue
            seen.add(S); out.append((c, S, self.kids(S)))
        return out
    def compat(self, A, patch):
        occ = set()
        for t in patch: occ |= self.cubes(t)
        pset = set(patch); out = []
        for (c, S, K) in self.supertiles_containing(A):
            ok = True
            for kd in K:
                if kd in pset: continue
                if not self.cubes(kd).isdisjoint(occ): ok = False; break
                if not all(self.pair_ok(kd, t) for t in patch): ok = False; break
            if ok: out.append((c, S, K))
        return out

    # ---- enclosures of `centre` compatible with fixed tiles
    def enclosures_given(self, centre, fixed):
        fixed = [t for t in fixed if t != centre]
        cand = [B for B in (self.absolute(centre, rp) for rp in self.rel_by_i[centre[2]])
                if B in fixed or all(self.pair_ok(B, f) for f in fixed)]
        M = len(cand)
        comp = [[True] * M for _ in range(M)]
        for i in range(M):
            for j in range(i + 1, M): comp[i][j] = comp[j][i] = self.pair_ok(cand[i], cand[j])
        cc = self.cubes(centre)
        panels = [(c, d) for c in cc for d in self.dirs if vadd(c, d) not in cc]
        covers = [frozenset(k for k, (c, d) in enumerate(panels) if vadd(c, d) in self.cubes(cand[i])) for i in range(M)]
        by_panel = [[i for i in range(M) if k in covers[i]] for k in range(len(panels))]
        out = []
        def bt(chosen, covered):
            k = 0
            while k < len(panels) and k in covered: k += 1
            if k == len(panels): out.append(tuple(cand[i] for i in chosen)); return
            for i in by_panel[k]:
                if all(comp[i][j] for j in chosen): bt(chosen + [i], covered | covers[i])
        bt([], frozenset()); return out

    # ---- liveness with parent propagation
    def enclosed(self, A, occ):
        ca = self.cubes(A)
        return all((vadd(c, d) in ca or vadd(c, d) in occ) for c in ca for d in self.dirs)
    def propagate(self, tiles, occ, decoded):
        changed = True
        while changed:
            changed = False
            for A in list(tiles):
                if A in decoded or not self.enclosed(A, occ): continue
                opts = self.compat(A, tiles)
                if not opts: return False
                if len(opts) == 1:
                    c, S, K = opts[0]; decoded[A] = S; pset = set(tiles)
                    for kd in K:
                        if kd not in pset:
                            tiles.append(kd); pset.add(kd)
                            for x in self.cubes(kd): occ[x] = kd
                            changed = True
        return True
    def live(self, patch, node_limit=50000):
        tiles = list(patch); occ = {}
        for T in tiles:
            for c in self.cubes(T): occ[c] = T
        decoded = {}
        if not self.propagate(tiles, occ, decoded): return False
        self.nodes = 0; self.limit_hit = False; self.node_limit = node_limit
        res = self._dfs(tiles, occ, decoded, list(patch))
        if self.limit_hit: self.limits = getattr(self, 'limits', 0) + 1; return True
        return res
    def _dfs(self, tiles, occ, decoded, must):
        self.nodes += 1
        if self.nodes > self.node_limit: self.limit_hit = True; return False
        best = None; bestc = None
        for A in must:
            ca = self.cubes(A)
            for c in ca:
                for d in self.dirs:
                    x = vadd(c, d)
                    if x in ca or x in occ: continue
                    cands = [B for B in self.candidates_for(A, c, d) if all(self.pair_ok(B, t) for t in tiles)]
                    if not cands: return False
                    if best is None or len(cands) < len(bestc):
                        best, bestc = (A, c, d), cands
                        if len(cands) == 1: break
                if best is not None and len(bestc) == 1: break
            if best is not None and len(bestc) == 1: break
        if best is None: return True
        for B in bestc:
            t2 = list(tiles); o2 = dict(occ); d2 = dict(decoded); t2.append(B)
            for x in self.cubes(B): o2[x] = B
            if not self.propagate(t2, o2, d2): continue
            if self._dfs(t2, o2, d2, must): return True
            if self.limit_hit: return False
        return False

# ---------------------------------------------------------------- languages and level rules
def contacts(MR, tiles):
    owner = {}
    for idx, T in enumerate(tiles):
        for c in MR.cubes(T): owner[c] = idx
    out = set()
    for idx, T in enumerate(tiles):
        for c in MR.cubes(T):
            for d in MR.dirs:
                j = owner.get(vadd(c, d))
                if j is not None and j != idx:
                    g, t = MR.rel(T, tiles[j]); out.add((g, t, T[2], tiles[j][2]))
    return out

def language_marked(ch):
    """R_1: contacts of marked tiles inside arbitrarily large supertiles (stabilised as in cert3.language)."""
    MR = MRules(ch, set(), True)
    S0 = (ch.id, tuple([0] * ch.N))
    K0 = MR.kids(S0)
    P1 = contacts(MR, K0)
    cache = {}
    def X(p):          # contacts between kids of S0 and kids of the neighbouring supertile 2p
        key = (p[0], p[1])
        if key in cache: return cache[key]
        g, t = key; S = (g, vscale(2, t))
        r = contacts(MR, K0 + MR.kids(S)); cache[key] = r; return r
    A = set(P1); lev = 1
    while True:
        B = set(P1)
        for p in A: B |= X(p)
        lev += 1
        if B == A: return frozenset(A), lev - 1
        A = B

def next_level_rule(MR):
    """R_{j+1}: relative poses (G,T) of supertiles pairwise admissible with S0 under MR (children rule)."""
    ch = MR.ch
    S0 = (ch.id, tuple([0] * ch.N)); K0 = MR.kids(S0)
    seen = set(); cand = set()
    for a in K0:
        for rp in MR.rel_by_i[a[2]]:
            b = MR.absolute(a, rp)
            for c in (range(MR.k)):           # b could be child c of its parent (any c if unmarked; must be b's index if marked)
                if MR.marked and c != b[2]: continue
                G0, T0 = MR.sup_for_child[c]
                S = (compose(b[0], G0), vadd(b[1], apply(b[0], T0)))
                if S in seen: continue
                seen.add(S)
                K = MR.kids(S); ok = True
                for x in K0:
                    for y in K:
                        if not MR.pair_ok(x, y): ok = False; break
                    if not ok: break
                if ok: cand.add(S)
    return cand

def claim1(MR, quiet=True, node_limit=50000):
    T = MR.origin; MR.limits = 0
    encl = MR.enclosures_given(T, [])
    st = dict(enclosures=len(encl), live=0, dead=0, U_fail=0, Pr_fail=0, typeA=0, typeB=0)
    for n, E in enumerate(encl):
        patch = [T] + list(E)
        if not MR.live(patch, node_limit): st['dead'] += 1; continue
        st['live'] += 1
        cs = MR.compat(T, patch)
        if len(cs) != 1: st['U_fail'] += 1; continue
        c, S, K_all = cs[0]
        missing = [K for K in K_all if K not in patch]
        if not missing: st['typeA'] += 1; continue
        st['typeB'] += 1
        ok_any = False
        for U in [K for K in K_all if K in patch and K != T]:
            good = True
            for EU in MR.enclosures_given(U, patch):
                joint = list(dict.fromkeys(patch + list(EU)))
                if all(K in joint for K in missing): continue
                if MR.live(joint, node_limit): good = False; break
            if good: ok_any = True; break
        if not ok_any: st['Pr_fail'] += 1
    st['limit'] = MR.limits
    return st

def run_marked(cells, dissection, tag, maxlevels=5, quiet=False, do_claims=True):
    t0 = time.time()
    cells = normalize_cells(cells)
    ch1 = level_tile(cells, dissection, 1)
    R1, lev = language_marked(ch1)
    out = dict(tag=tag, R1=len(R1), lang_level=lev, levels=[])
    if not quiet: log(f"=== {tag}: |R1|={len(R1)} (marked, stabilised at level {lev}); unmarked projection {len({(g,t) for g,t,_,_ in R1})}")
    rules = [None, R1]; MRs = [None, MRules(ch1, R1, True)]
    status = 'open'
    for j in range(1, maxlevels + 1):
        MR = MRs[j]
        cand = next_level_rule(MR)
        odd = [q for q in cand if any(x % (2 ** j) for x in q[1])]
        if odd:
            status = f'offset@{j}({len(odd)})'; out['levels'].append((j, len(cand), 'offset'))
            if not quiet: log(f"  level {j}->{j+1}: {len(cand)} supertile poses, {len(odd)} OFFSET -> fail")
            break
        Rn = frozenset((g, t, 0, 0) for (g, t) in cand)   # unit cells
        prev_unmarked = frozenset((g, vscale(2, t), 0, 0) for (g, t, _, _) in rules[j])
        same = (Rn == prev_unmarked)
        if not quiet: log(f"  level {j}->{j+1}: |R_{j+1}|={len(Rn)} (prev unmarked {len(prev_unmarked)}) self-similar={same} ({time.time()-t0:.0f}s)")
        out['levels'].append((j, len(Rn), same))
        if same and j >= 2:
            status = f'self-similar@{j}'; break
        rules.append(Rn); chj = level_tile(cells, dissection, j + 1); MRs.append(MRules(chj, Rn, False))
    out['status'] = status
    if not status.startswith('self-similar') or not do_claims:
        if not quiet: log(f"=== {tag}: {status} ({time.time()-t0:.0f}s)")
        return out
    # enclosure analysis at every level 1..j
    allok = True
    for jj in range(1, len(MRs)):
        st = claim1(MRs[jj], quiet=quiet)
        out[f'claim1_L{jj}'] = st
        ok = st['U_fail'] == 0 and st['Pr_fail'] == 0 and st['limit'] == 0
        if not quiet: log(f"  claim1 level {jj}: {st} ok={ok} ({time.time()-t0:.0f}s)")
        allok &= ok
    out['certified'] = allok
    if not quiet: log(f"=== {tag}: {status} CERTIFIED={allok} ({time.time()-t0:.0f}s)")
    return out

if __name__ == '__main__':
    # validation: the chair with Paper-1 frames
    from chairN import Chair
    fr30 = {'0': (0, 1, 2), (0, 0, 0): (0, 1, 2), (0, 0, 1): (0, 2, 1), (0, 1, 0): (2, 1, 0), (0, 1, 1): (1, 2, 0),
            (1, 0, 0): (1, 0, 2), (1, 0, 1): (0, 1, 2), (1, 1, 0): (2, 0, 1)}
    c = Chair(3, fr30); dis = []
    for (h, tau) in c.children:
        s = apply(h, (1, 1, 1)); dis.append((h, tuple(tau[i] - (1 - s[i]) // 2 for i in range(3))))
    cells = frozenset(u for u in itertools.product((0, 1), repeat=3) if u != (1, 1, 1))
    r = run_marked(cells, dis, 'chair_marked')
    print(r['status'], r.get('certified'))
