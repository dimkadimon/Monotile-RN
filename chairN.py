"""
chairN.py -- Finite certificate machinery for the N-dimensional chair (corner-deleted hypercube)
with full-frame facet matching rules ("maximal rules" in the sense of Tsiokos / Goodman-Strauss).

Conventions
-----------
* A unit cube is identified by its centre, an odd integer vector in Z^N.
* C_N = [0,2]^N minus (1,2]^N has cube centres {2u+1 : u in {0,1}^N, u != (1,..,1)}.
* An element g of the hyperoctahedral group B_N is a signed permutation, stored as
  (perm, signs) with g(x)_i = signs[i] * x[perm[i]].
* A tile pose is (g, t) with t in Z^N: cubes = { g(2u+1) + 2t }.
* Level-1 supertile in the local frame: 2C_N = [0,4]^N \ (2,4]^N, children
      T_0 = (pi_0, (1,..,1)),  T_v = (D_s pi_v, 4v)  for v in {0,1}^N \ {1},  s_i = +1 if v_i = 0 else -1,
  where pi_0, pi_v are coordinate permutations (symmetries of C_N) -- the *design freedom*.
"""
import itertools, sys
from collections import defaultdict

# ---------------------------------------------------------------- group B_N
def make_group(N):
    els = []
    for perm in itertools.permutations(range(N)):
        for signs in itertools.product((1, -1), repeat=N):
            els.append((perm, signs))
    return els

def apply(g, x):
    perm, signs = g
    return tuple(signs[i] * x[perm[i]] for i in range(len(perm)))

def compose(g, h):            # (g*h)(x) = g(h(x))
    pg, sg = g; ph, sh = h
    N = len(pg)
    # h(x)_j = sh[j] x[ph[j]] ; g(y)_i = sg[i] y[pg[i]] = sg[i] sh[pg[i]] x[ph[pg[i]]]
    return (tuple(ph[pg[i]] for i in range(N)), tuple(sg[i] * sh[pg[i]] for i in range(N)))

def inverse(g):
    perm, signs = g
    N = len(perm)
    inv = [0] * N; isg = [1] * N
    for i in range(N):
        inv[perm[i]] = i
        isg[perm[i]] = signs[i]
    return (tuple(inv), tuple(isg))

def det(g):
    perm, signs = g
    N = len(perm)
    # sign of permutation
    seen = [False] * N; s = 1
    for i in range(N):
        if not seen[i]:
            j = i; L = 0
            while not seen[j]:
                seen[j] = True; j = perm[j]; L += 1
            if L % 2 == 0: s = -s
    for x in signs: s *= x
    return s

def identity(N):
    return (tuple(range(N)), tuple([1] * N))

def perm_only(perm):
    return (tuple(perm), tuple([1] * len(perm)))

def diag(signs):
    return (tuple(range(len(signs))), tuple(signs))

def vadd(a, b): return tuple(x + y for x, y in zip(a, b))
def vsub(a, b): return tuple(x - y for x, y in zip(a, b))
def vscale(k, a): return tuple(k * x for x in a)

# ---------------------------------------------------------------- the tile
class Chair:
    def __init__(self, N, frames=None):
        """frames: dict mapping child key ('0' or v-tuple) -> coordinate permutation (tuple)."""
        self.N = N
        self.G = make_group(N)
        self.id = identity(N)
        self.U = [u for u in itertools.product((0, 1), repeat=N) if any(c == 0 for c in u)]
        self.Uset = set(self.U)
        self.local_cubes = [tuple(2 * c + 1 for c in u) for u in self.U]
        self.dirs = []
        for i in range(N):
            for s in (1, -1):
                d = [0] * N; d[i] = 2 * s; self.dirs.append(tuple(d))
        # local panels: (cube centre, direction) with neighbour cube outside the tile
        cs = set(self.local_cubes)
        self.panels = [(c, d) for c in self.local_cubes for d in self.dirs if vadd(c, d) not in cs]
        self.panel_index = {p: i for i, p in enumerate(self.panels)}
        assert len(self.panels) == N * 2 ** N
        # children of the level-1 supertile
        if frames is None:
            frames = {}
        ident = tuple(range(N))
        self.children = []  # list of (h, tau)
        pi0 = perm_only(frames.get('0', ident))
        self.children.append((pi0, tuple([1] * N)))
        for v in self.U:
            s = tuple(1 if c == 0 else -1 for c in v)
            h = compose(diag(s), perm_only(frames.get(v, ident)))
            self.children.append((h, tuple(4 * c for c in v)))
        # sanity: children partition 2C_N
        allc = []
        for (h, tau) in self.children:
            allc += self.cubes((h, tau))
        big = {tuple(2 * c + 1 for c in u) for u in itertools.product((0, 1, 2, 3), repeat=N)
               if any(c < 2 for c in u)}
        assert len(allc) == len(set(allc)) == len(big) and set(allc) == big, "children do not partition 2C_N"

    def cubes(self, pose):
        g, t = pose
        tt = vscale(2, t)
        return [vadd(apply(g, c), tt) for c in self.local_cubes]

    def sub_children(self, G, T, k):
        """children (level k-1) of a level-k supertile with pose (G,T)."""
        out = []
        for (h, tau) in self.children:
            out.append((compose(G, h), vadd(T, vscale(2 ** (k - 1), apply(G, tau)))))
        return out

    def patch(self, k):
        """all tiles of the level-k supertile with pose (id,0)."""
        level = [(self.id, tuple([0] * self.N))]
        for lev in range(k, 0, -1):
            nxt = []
            for (G, T) in level:
                nxt += self.sub_children(G, T, lev)
            level = nxt
        return level

    # relative pose of B in the frame of A
    def rel(self, A, B):
        gA, tA = A; gB, tB = B
        gi = inverse(gA)
        return (compose(gi, gB), apply(gi, vsub(tB, tA)))

    def contacts(self, tiles):
        """ordered relative poses (B rel A) for all panel-adjacent tile pairs in a list of tiles."""
        owner = {}
        for idx, T in enumerate(tiles):
            for c in self.cubes(T):
                owner[c] = idx
        poses = set()
        for idx, T in enumerate(tiles):
            for c in self.cubes(T):
                for d in self.dirs:
                    j = owner.get(vadd(c, d))
                    if j is not None and j != idx:
                        poses.add(self.rel(T, tiles[j]))
        return poses

    def panel_triples(self, pose):
        """facet-level contacts between tile (id,0) and tile at `pose`:
        set of (panelA_index, panelB_index, g_rel)."""
        g, t = pose
        gi = inverse(g)
        Bc = set(self.cubes(pose))
        trip = set()
        for (c, d) in self.panels:
            nb = vadd(c, d)
            if nb in Bc:
                # local coords of nb in B: gi(nb - 2t); local direction of the panel of B: -gi(d)
                lc = apply(gi, vsub(nb, vscale(2, t)))
                ld = apply(gi, tuple(-x for x in d))
                pb = self.panel_index[(lc, ld)]
                trip.add((self.panel_index[(c, d)], pb, g))
        return trip

    def overlaps(self, pose):
        return not set(self.cubes(pose)).isdisjoint(self.local_cubes)

    def t_ranges(self, img, hi_target):
        """per-coordinate ranges of t such that the shifted centres [lo+2t, hi+2t] can touch [-1, hi_target+... ]"""
        rs = []
        for i in range(self.N):
            lo = min(c[i] for c in img); hi = max(c[i] for c in img)
            # target centres occupy [1, hi_target]; touching requires overlap of [lo+2t,hi+2t] with [-1, hi_target+2]
            tmin = -((hi + 1) // 2) - 1; tmax = (hi_target + 2 - lo) // 2 + 1
            rs.append(range(tmin, tmax + 1))
        return rs

    def all_touching_poses(self):
        """all poses (g,t) with tile disjoint from (id,0) and sharing >= 1 panel with it."""
        base = set(self.local_cubes)
        out = []
        for g in self.G:
            img = [apply(g, c) for c in self.local_cubes]
            for t in itertools.product(*self.t_ranges(img, 3)):
                tt = vscale(2, t)
                cs = [vadd(c, tt) for c in img]
                if any(c in base for c in cs):
                    continue
                cset = set(cs)
                touch = False
                for c in self.local_cubes:
                    for d in self.dirs:
                        if vadd(c, d) in cset:
                            touch = True; break
                    if touch: break
                if touch:
                    out.append((g, t))
        return out
