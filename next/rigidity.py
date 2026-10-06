"""
rigidity.py -- "mixed substitution" construction: the rigidity theorem for ancestry markings.

Setting.  T polycube, D a 2-dissection (2T = union of 2^N copies), hierarchy = D at every level.
Depth-k ancestry marking: a tile carries its frame and the indices (i_1..i_k) of its ancestors
(i_j = index of its level-j ancestor inside its level-(j+1) ancestor).  L_k = set of contacts
(g_rel, t_rel, mark_A, mark_B) between face-adjacent marked tiles in the D-hierarchy.  k=0: frames only.

Mixed tiling  T*(D', M):  a D-hierarchical tiling in which every level-M supertile (shape 2^(M-1) T) is
dissected by D' instead of D (its level-(M-1) children are pure D).  For M >= k+2 all marks are honest.

Theorem (mixed substitution).  If
  (a) every contact inside a D'-dissected level-M supertile lies in L_k, and
  (b) for every relative pose p in the adjacency language P of the D-hierarchy, every contact between
      two D'-dissected level-M supertiles at relative pose 2^(M-1) p lies in L_k,
then T*(D', M) is a tiling admitted by the depth-k contact rule.  If moreover (k >= 1, M = k+2)
  (c) D' places two of its children at a relative pose not in P,   or
  (c') some mark-group (level-(k+1) supertile, readable from the marks) of T* has no D-parent all of
       whose children are mark-groups of T*,
then T* is not in the hull of the marked D-hierarchy: the rule does not enforce the hierarchy.

Proof sketch.  (a)+(b) cover all adjacent pairs: pairs inside one level-M supertile, and pairs in two
different level-M supertiles, whose relative pose is in 2^(M-1) P because the hierarchy above level M
is untouched.  In a marked D-hierarchical tiling the marks are honest, so mark-groups are true
level-(k+1) supertiles, adjacent groups have poses in P, and every group has a D-parent whose children
are groups; (c)/(c') violate this.                                                      QED
"""
import sys, time, itertools, pickle
sys.path.insert(0, '/home/user/research'); sys.path.insert(0, '/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from chairN import apply, compose, inverse, vadd, vsub, vscale, identity, make_group
from polytile import PolyTile, normalize_cells, scaled_cells, dissections, dissection_classes, stabiliser, frame_variants
from mcert import level_tile
from cert3 import log

class Hier:
    def __init__(self, cells, D):
        self.cells = normalize_cells(cells); self.N = len(next(iter(self.cells)))
        self.D = list(D); self.ch1 = PolyTile(self.cells, self.D)
        self.id = identity(self.N); self.dirs = self.ch1.dirs
        self._kids = {}
        self._cubes = {}
    def kids_maps(self, dis, m):
        """(h, tau) child maps of a level-m supertile for dissection dis (level-(m-1) children)."""
        key = (tuple(dis), m)
        if key not in self._kids: self._kids[key] = level_tile(self.cells, dis, m - 1).children
        return self._kids[key]
    def cubes(self, g, t):
        key = (g, t); r = self._cubes.get(key)
        if r is None:
            r = frozenset(vadd(apply(g, c), vscale(2, t)) for c in self.ch1.local_cubes); self._cubes[key] = r
        return r
    def expand(self, S, m, k, dis_at, path=()):
        """tiles (g, t, mark) of level-m supertile S=(G,T); dis_at(level) -> dissection used at that level.
        mark = indices of ancestors at levels 2..k+1 (needs m >= k+1)."""
        G, T = S
        if m == 1:
            mk = tuple(reversed(path))[:k]      # path: indices at levels m..2 (top first)
            assert len(mk) == k
            return [(G, T, mk)]
        out = []
        for i, (h, tau) in enumerate(self.kids_maps(dis_at(m), m)):
            out += self.expand((compose(G, h), vadd(T, apply(G, tau))), m - 1, k, dis_at, path + (i,))
        return out
    def contacts(self, tiles):
        owner = {}
        for idx, (g, t, mk) in enumerate(tiles):
            for c in self.cubes(g, t): owner[c] = idx
        out = set()
        for idx, (g, t, mk) in enumerate(tiles):
            gi = inverse(g)
            for c in self.cubes(g, t):
                for d in self.dirs:
                    j = owner.get(vadd(c, d))
                    if j is not None and j != idx:
                        g2, t2, mk2 = tiles[j]
                        out.add((compose(gi, g2), apply(gi, vsub(t2, t)), mk, mk2))
        return out
    # adjacency language of the D-hierarchy (frame-marked, unmarked indices)
    def adjacency_language(self):
        pureD = lambda m: self.D
        S0 = (self.id, tuple([0] * self.N))
        K0 = self.expand(S0, 2, 0, pureD)
        P1 = {(g, t) for (g, t, _, _) in self.contacts(K0)}
        cache = {}
        def X(p):
            if p in cache: return cache[p]
            g, t = p
            r = {(a, b) for (a, b, _, _) in self.contacts(K0 + self.expand((g, vscale(2, t)), 2, 0, pureD))}
            cache[p] = r; return r
        A = set(P1); lev = 1
        while True:
            B = set(P1)
            for p in A: B |= X(p)
            lev += 1
            if B == A: self.P = frozenset(A); self.Plev = lev - 1; return self.P
            A = B
    def language(self, k):
        """L_k: contacts of depth-k marked tiles in the D-hierarchy."""
        if not hasattr(self, 'P'): self.adjacency_language()
        pureD = lambda m: self.D
        S0 = (self.id, tuple([0] * self.N)); m = k + 1
        K0 = self.expand(S0, m, k, pureD)
        L = set(self.contacts(K0))
        f = 2 ** (m - 1)
        for (g, t) in self.P:
            L |= self.contacts(K0 + self.expand((g, vscale(f, t)), m, k, pureD))
        return frozenset(L)
    def mixed_check(self, Dp, k, M, L=None):
        """conditions (a), (b) for T*(D', M) under L_k; returns (ok, n_bad_inside, n_bad_between, groups info)."""
        if L is None: L = self.language(k)
        if not hasattr(self, 'P'): self.adjacency_language()
        dis_at = lambda m: (list(Dp) if m == M else self.D)
        S0 = (self.id, tuple([0] * self.N))
        K0 = self.expand(S0, M, k, dis_at)
        bad_in = {c for c in self.contacts(K0) if c not in L}
        bad_between = 0
        f = 2 ** (M - 1)
        if not bad_in:
            for (g, t) in self.P:
                K1 = self.expand((g, vscale(f, t)), M, k, dis_at)
                bad_between += sum(1 for c in self.contacts(K0 + K1) if c not in L)
        return (not bad_in and bad_between == 0), len(bad_in), bad_between
    def nonhier_c(self, Dp):
        """criterion (c): D' child pairs (unscaled poses, as level-1 copies) adjacent with pose not in P."""
        if not hasattr(self, 'P'): self.adjacency_language()
        S0 = (self.id, tuple([0] * self.N))
        K = self.expand(S0, 2, 0, lambda m: list(Dp))
        pairs = {(g, t) for (g, t, _, _) in self.contacts(K)}
        return {p for p in pairs if p not in self.P}

def same_partition(cells, D1, D2):
    f = lambda D: frozenset(frozenset(vadd(apply(g, c), t) for c in cells) for (g, t) in D)
    return f(D1) == f(D2)

def analyse(name, cells, D, k, M, Dprimes, quiet=False):
    H = Hier(cells, D); P = H.adjacency_language(); L = H.language(k)
    if not quiet: log(f"[{name}] |P|={len(P)} (stable at level {H.Plev}); |L_{k}|={len(L)}; testing {len(Dprimes)} D' at M={M}")
    res = []
    for di, Dp in enumerate(Dprimes):
        ok, bi, bb = H.mixed_check(Dp, k, M, L)
        c = H.nonhier_c(Dp) if ok else set()
        res.append(dict(di=di, ok=ok, bad_in=bi, bad_between=bb, nonhier_pairs=len(c), same=same_partition(H.cells, D, Dp)))
        if not quiet: log(f"  D'#{di}: admissible={ok} (bad inside {bi}, between {bb}); same partition as D: {res[-1]['same']}; (c) pairs not in P: {len(c)}")
    return H, res

if __name__ == '__main__':
    from census3 import SHAPES
    name = sys.argv[1]; k = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    M = int(sys.argv[3]) if len(sys.argv) > 3 else k + 2
    maxc = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 9
    cells = normalize_cells(frozenset(SHAPES[name])); G = make_group(3)
    classes, ntot, nsym = dissection_classes(cells, G)
    allD = dissections(cells, G)
    log(f"### {name}: {ntot} dissections, {len(classes)} classes; depth k={k}, D' level M={M}")
    summary = []
    for ci, D in enumerate(classes[:maxc]):
        Dprimes = [Dp for Dp in allD if not same_partition(cells, D, Dp)]
        # also frame variants of D itself (same partition, different child frames), excluding D
        stab = stabiliser(cells, G); fv = []
        for i in range(1, len(D)):
            for (s_, u) in stab[1:]:
                g, t = D[i]; Dp = list(D); Dp[i] = (compose(g, s_), vadd(t, apply(g, u))); fv.append(Dp)
        H, res = analyse(f"{name} class {ci}", cells, D, k, M, Dprimes + fv, quiet=True)
        n_ok = sum(r['ok'] for r in res); n_nh = sum(r['ok'] and r['nonhier_pairs'] > 0 for r in res)
        log(f"class {ci}: |P|={len(H.P)} |L_{k}|={len(H.language(k))}; D' candidates {len(res)} (other partitions {len(Dprimes)}, frame variants {len(fv)}): admissible {n_ok}, admissible & (c) {n_nh}")
        summary.append((ci, D, res))
        pickle.dump(summary, open(f'rigidity_{name}_k{k}_M{M}.pkl', 'wb'))
