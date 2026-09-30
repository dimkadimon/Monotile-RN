"""
cert3.py -- optimised, generic finite certificate.

Generic "shape + allowed relative poses" machinery (class Rules) used twice:
  * tile level:      shape = C_N,  allowed = P  (contact language, tight w.r.t. facet rules F)
  * supertile level: shape = 2C_N, allowed = P_S (pairwise child-admissible supertile poses)

Claim 1 (tile level):  every 2-extendable complete enclosure of a tile is type A or type B (see certify2.py).
Claim 2 (supertile level): every 2-extendable complete enclosure of a supertile by pairwise-admissible
                            supertiles uses only poses in 2P.
Stabilisation lemma: A_{k+1} = P_1 ∪ X(A_k) is monotone, so A_k = A_{k+1} proves the language is stable forever.
"""
import itertools, time, sys, pickle, os
from chairN import *

def log(*a):
    print(*a, flush=True)

class Rules:
    def __init__(self, ch, local_cubes, P):
        self.ch = ch; self.N = ch.N
        self.local = tuple(local_cubes)
        self.local_set = frozenset(local_cubes)
        self.P = frozenset(P); self.Plist = sorted(P)
        self.dirs = ch.dirs
        self._cubes = {}; self._halo = {}
        # local panels & which rel poses cover each
        self.panels = [(c, d) for c in self.local for d in self.dirs if vadd(c, d) not in self.local_set]
        self.pidx = {p: i for i, p in enumerate(self.panels)}
        self.cover = [[] for _ in self.panels]
        self.covers_of = {}
        origin = (ch.id, tuple([0] * self.N))
        for rp in self.Plist:
            cb = self.cubes(rp)
            ks = frozenset(i for i, (c, d) in enumerate(self.panels) if vadd(c, d) in cb)
            self.covers_of[rp] = ks
            for k in ks: self.cover[k].append(rp)

    def cubes(self, pose):
        r = self._cubes.get(pose)
        if r is None:
            g, t = pose; tt = vscale(2, t)
            r = frozenset(vadd(apply(g, c), tt) for c in self.local)
            self._cubes[pose] = r
        return r

    def halo(self, pose):
        r = self._halo.get(pose)
        if r is None:
            cb = self.cubes(pose)
            r = frozenset(vadd(c, d) for c in cb for d in self.dirs) - cb
            self._halo[pose] = r
        return r

    def absolute(self, A, rp):
        gA, tA = A; g, t = rp
        return (compose(gA, g), vadd(tA, apply(gA, t)))

    def pair_ok(self, A, B):
        ca = self.cubes(A); cb = self.cubes(B)
        if not ca.isdisjoint(cb): return False
        if self.halo(A).isdisjoint(cb): return True
        return self.ch.rel(A, B) in self.P

    def local_panel_index(self, A, c, d):
        gA, tA = A; gi = inverse(gA)
        lc = apply(gi, vsub(c, vscale(2, tA))); ld = apply(gi, d)
        return self.pidx[(lc, ld)]

    def candidates_for(self, A, c, d):
        k = self.local_panel_index(A, c, d)
        return [self.absolute(A, rp) for rp in self.cover[k]]

    def open_panels(self, A, occ):
        ca = self.cubes(A)
        return [(A, c, d) for c in ca for d in self.dirs if vadd(c, d) not in ca and vadd(c, d) not in occ]

    # ------------------------------------------------------------ enclosures of a single shape
    def enclosures(self, centre):
        cand = [self.absolute(centre, rp) for rp in self.Plist]
        M = len(cand)
        comp = [[True] * M for _ in range(M)]
        for i in range(M):
            for j in range(i + 1, M):
                comp[i][j] = comp[j][i] = self.pair_ok(cand[i], cand[j])
        cc = self.cubes(centre)
        panels = [(c, d) for c in cc for d in self.dirs if vadd(c, d) not in cc]
        covers = [frozenset(k for k, (c, d) in enumerate(panels) if vadd(c, d) in self.cubes(cand[i])) for i in range(M)]
        by_panel = [[i for i in range(M) if k in covers[i]] for k in range(len(panels))]
        out = []
        def bt(chosen, covered):
            k = 0
            while k < len(panels) and k in covered: k += 1
            if k == len(panels):
                out.append(tuple(cand[i] for i in chosen)); return
            for i in by_panel[k]:
                if all(comp[i][j] for j in chosen):
                    bt(chosen + [i], covered | covers[i])
        bt([], frozenset())
        return out

    # ------------------------------------------------------------ extendability
    def lookahead_ok(self, tiles, occ, must):
        """1-step: every open panel of tiles in `must` has at least one consistent candidate."""
        for A in must:
            for (_, c, d) in self.open_panels(A, occ):
                if not any(all(self.pair_ok(B, T) for T in tiles) for B in self.candidates_for(A, c, d)):
                    return False
        return True

    def extend_exists(self, tiles, must, node_limit=100000):
        occ = {}
        for T in tiles:
            for c in self.cubes(T): occ[c] = T
        need = []
        for A in must:
            need += self.open_panels(A, occ)
        self.nodes = 0; self.limit_hit = False
        return self._dfs(list(tiles), occ, need, node_limit)

    def _dfs(self, tiles, occ, need, node_limit):
        self.nodes += 1
        if self.nodes > node_limit:
            self.limit_hit = True; return False
        # choose the uncovered panel with fewest consistent candidates (fail-first)
        best = None; bestc = None
        for (A, c, d) in need:
            if vadd(c, d) in occ: continue
            cands = [B for B in self.candidates_for(A, c, d) if all(self.pair_ok(B, T) for T in tiles)]
            if not cands: return False
            if best is None or len(cands) < len(bestc):
                best, bestc = (A, c, d), cands
                if len(cands) == 1: break
        if best is None: return True
        for B in bestc:
            cb = self.cubes(B)
            for x in cb: occ[x] = B
            tiles.append(B)
            if self._dfs(tiles, occ, need, node_limit):
                tiles.pop()
                for x in cb: del occ[x]
                return True
            tiles.pop()
            for x in cb: del occ[x]
        return False


# ---------------------------------------------------------------- language
def language(ch, verbose=True):
    """A_1 = P_1 ; A_{k+1} = P_1 ∪ X(A_k). Returns stable set and the level at which it stabilised."""
    P1 = ch.contacts(ch.patch(1))
    S0 = ch.sub_children(ch.id, tuple([0] * ch.N), 1)
    cache = {}
    def X(p):
        if p in cache: return cache[p]
        g, t = p
        S1 = [(compose(g, h), vadd(vscale(2, t), apply(g, tau))) for (h, tau) in ch.children]
        r = ch.contacts(S0 + S1)
        cache[p] = r; return r
    A = set(P1); k = 1
    while True:
        B = set(P1)
        for p in A: B |= X(p)
        k += 1
        if verbose: log(f"  level {k}: |A|={len(B)}")
        if B == A: return frozenset(A), k - 1
        A = B

def facet_rules(ch, P):
    F = set()
    for p in P: F |= ch.panel_triples(p)
    return F

def tightness(ch, P, F):
    allp = ch.all_touching_poses()
    PF = set(p for p in allp if (lambda tr: tr and tr <= F)(ch.panel_triples(p)))
    return len(allp), PF

# ---------------------------------------------------------------- claims
def claim1(ch, R, node_limit=100000):
    N = ch.N
    centre = (ch.id, tuple([0] * N))
    children = ch.children
    def sup_for_child(c):
        h, tau = children[c]; G = inverse(h)
        return (G, tuple(-x for x in apply(G, tau)))
    def kids(G, T):
        return [(compose(G, h), vadd(T, apply(G, tau))) for (h, tau) in children]
    central_rel = {}
    for c in range(1, len(children)):
        central_rel[kids(*sup_for_child(c))[0]] = c
    central_kids = set(kids(*sup_for_child(0))) - {centre}
    t0 = time.time()
    encl = R.enclosures(centre)
    log(f"  Claim1: complete enclosures = {len(encl)}  ({time.time()-t0:.0f}s)")
    st = dict(total=len(encl), typeA=0, typeB=0, dead1=0, dead2=0, bad=0, c4_fail=0, limit=0)
    bad = []
    for n, E in enumerate(encl):
        tiles = set(E)
        isA = central_kids <= tiles
        Ecubes = set().union(*[R.cubes(t) for t in E]) | R.cubes(centre)
        compatU = []
        for u in tiles:
            if u not in central_rel: continue
            ok = True
            for kd in kids(*sup_for_child(central_rel[u])):
                if kd == centre or kd in tiles: continue
                if not R.cubes(kd).isdisjoint(Ecubes): ok = False; break
            if ok: compatU.append(u)
        good = (isA and not compatU) or ((not isA) and len(compatU) == 1)
        c4ok = (not any(u in central_kids for u in tiles)) or isA
        if good and c4ok:
            st['typeA' if isA else 'typeB'] += 1; continue
        all_t = [centre] + list(E)
        occ = {}
        for T in all_t:
            for c in R.cubes(T): occ[c] = T
        if not R.lookahead_ok(all_t, occ, all_t):
            st['dead1'] += 1; continue
        ext = R.extend_exists(all_t, all_t, node_limit)
        if R.limit_hit: st['limit'] += 1
        if ext:
            st['bad'] += 1
            if not c4ok: st['c4_fail'] += 1
            bad.append(E)
        else:
            st['dead2'] += 1
        if n % 100 == 0: log(f"    ...{n}/{len(encl)} {st}")
    st['time'] = time.time() - t0
    log("  Claim1 stats:", st)
    return st, bad

def supertile_rules(ch, R_tile):
    """pairwise child-admissible supertile poses (relative to supertile at (id,0)), as a Rules object."""
    N = ch.N
    children = ch.children
    def kids(G, T):
        return [(compose(G, h), vadd(T, apply(G, tau))) for (h, tau) in children]
    S0 = kids(ch.id, tuple([0] * N))
    local = sorted(set().union(*[R_tile.cubes(k) for k in S0]))
    S0c = frozenset(local)
    S0h = frozenset(vadd(c, d) for c in S0c for d in ch.dirs) - S0c
    cand = set()
    for G in ch.G:
        img = [apply(G, c) for c in local]
        for T in itertools.product(*ch.t_ranges(img, 4 * 2 - 1)):
            tt = vscale(2, T)
            cs = frozenset(vadd(c, tt) for c in img)
            if not cs.isdisjoint(S0c) or S0h.isdisjoint(cs): continue
            # children-level check
            K = kids(G, T)
            ok = True
            for a in S0:
                for b in K:
                    if not R_tile.pair_ok(a, b): ok = False; break
                if not ok: break
            if ok: cand.add((G, T))
    return Rules(ch, local, cand), cand

def claim2(ch, R_tile, P, node_limit=100000):
    t0 = time.time()
    RS, cand = supertile_rules(ch, R_tile)
    scaledP = set((g, vscale(2, t)) for (g, t) in P)
    log(f"  Claim2: pairwise-admissible supertile poses = {len(cand)}; 2P={len(scaledP)}, 2P⊆cand: {scaledP <= cand}  ({time.time()-t0:.0f}s)")
    centre = (ch.id, tuple([0] * ch.N))
    encl = RS.enclosures(centre)
    log(f"  Claim2: complete supertile enclosures = {len(encl)}  ({time.time()-t0:.0f}s)")
    st = dict(total=len(encl), good=0, dead1=0, dead2=0, bad=0, limit=0)
    bad = []
    for n, E in enumerate(encl):
        if all(p in scaledP for p in E):
            st['good'] += 1; continue
        all_t = [centre] + list(E)
        occ = {}
        for T in all_t:
            for c in RS.cubes(T): occ[c] = T
        if not RS.lookahead_ok(all_t, occ, all_t):
            st['dead1'] += 1; continue
        ext = RS.extend_exists(all_t, all_t, node_limit)
        if RS.limit_hit: st['limit'] += 1
        if ext: st['bad'] += 1; bad.append(E)
        else: st['dead2'] += 1
        if n % 200 == 0: log(f"    ...{n}/{len(encl)} {st}")
    st['time'] = time.time() - t0
    log("  Claim2 stats:", st)
    return st, bad

def run(N, frames, tag):
    ch = Chair(N, frames)
    log(f"=== N={N} frames={frames}")
    t0 = time.time()
    P, k = language(ch)
    F = facet_rules(ch, P)
    ntouch, PF = tightness(ch, P, F)
    log(f"  |P|={len(P)} stable at level {k}; |F|={len(F)}; touching={ntouch}; |P_F|={len(PF)}; tight={PF == P}  ({time.time()-t0:.0f}s)")
    out = dict(N=N, frames=frames, P=sorted(P), F=sorted(F), stable_level=k, touching=ntouch, tight=(PF == P))
    if PF != P:
        out['PF_extra'] = sorted(PF - P)
        pickle.dump(out, open(f'cert_{tag}.pkl', 'wb')); return out
    R = Rules(ch, ch.local_cubes, P)
    s1, b1 = claim1(ch, R)
    out.update(claim1=s1, bad1=b1)
    pickle.dump(out, open(f'cert_{tag}.pkl', 'wb'))
    s2, b2 = claim2(ch, R, P)
    out.update(claim2=s2, bad2=b2)
    pickle.dump(out, open(f'cert_{tag}.pkl', 'wb'))
    log(f"=== done {tag}: claim1 bad={s1['bad']} c4_fail={s1['c4_fail']} limit={s1['limit']} | claim2 bad={s2['bad']} limit={s2['limit']}  total {time.time()-t0:.0f}s")
    return out

if __name__ == "__main__":
    N = int(sys.argv[1]); frames = eval(sys.argv[2]); tag = sys.argv[3]
    run(N, frames, tag)
