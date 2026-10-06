"""
polytile.py -- generalisation of chairN.Chair to an arbitrary polycube rep-2^N tile with a given dissection.

A PolyTile quacks like chairN.Chair (same attributes/methods used by cert3.Rules, cert3.language,
cert3.facet_rules and cert6.HCert / closure / claim1):
    N, G, id, local_cubes, dirs, panels, panel_index, children,
    cubes(pose), sub_children(G,T,k), patch(k), rel(A,B), contacts(tiles), panel_triples(pose),
    overlaps(pose), all_touching_poses()

Conventions (identical to chairN):
  * unit cube <-> its centre, an odd integer vector: cell x in Z^N  <->  centre 2x+1
  * pose (g, t), g in B_N (signed permutation), t in Z^N: cubes = g(centre) + 2t
  * level-1 supertile = 2*T in the local frame, children (h_i, tau_i) with cubes h_i(local) + 2 tau_i

A dissection is a list of (g_i, t_i): child i occupies g_i(cells) + t_i (cell coordinates inside the
scaled tile 2*cells); frames of symmetric tiles are chosen by right-multiplying g_i with stabiliser
elements (see `frame_variants`).
"""
import sys, itertools
sys.path.insert(0, '/home/user/research')
from chairN import make_group, apply, compose, inverse, det, identity, vadd, vsub, vscale

def normalize_cells(cells):
    N = len(next(iter(cells)))
    m = [min(c[i] for c in cells) for i in range(N)]
    return frozenset(tuple(c[i] - m[i] for i in range(N)) for c in cells)

def scaled_cells(cells, k=2):
    N = len(next(iter(cells)))
    return frozenset(tuple(k * c[i] + d[i] for i in range(N)) for c in cells
                     for d in itertools.product(range(k), repeat=N))

def stabiliser(cells, G):
    """elements g of G with normalize(g(cells)) == normalize(cells), returned with the translation
    s.t. g(cells)+t == cells  (as (g, t) in cell coordinates)."""
    cells = normalize_cells(cells); out = []
    for g in G:
        img = [apply(g, c) for c in cells]
        if normalize_cells(img) == cells:
            m = [min(c[i] for c in img) for i in range(len(g[0]))]
            out.append((g, tuple(-x for x in m)))
    return out

class PolyTile:
    def __init__(self, cells, dissection, name='tile', scale=2):
        cells = normalize_cells(cells)
        self.name = name; self.scale = scale
        self.cells = cells
        N = len(next(iter(cells))); self.N = N
        self.G = make_group(N); self.id = identity(N)
        self.local_cubes = [tuple(2 * x + 1 for x in c) for c in sorted(cells)]
        self.local_set = frozenset(self.local_cubes)
        self.dirs = []
        for i in range(N):
            for s in (1, -1):
                d = [0] * N; d[i] = 2 * s; self.dirs.append(tuple(d))
        self.panels = [(c, d) for c in self.local_cubes for d in self.dirs if vadd(c, d) not in self.local_set]
        self.panel_index = {p: i for i, p in enumerate(self.panels)}
        self.hi = [max(c[i] for c in self.local_cubes) for i in range(N)]
        # children: (h, tau) with h(2c+1) + 2 tau = 2(g c + t) + 1  =>  tau = t + (1 - g(1))/2
        self.dissection = list(dissection)
        self.children = []
        ones = tuple([1] * N)
        for (g, t) in dissection:
            s = apply(g, ones)                      # = signs vector
            tau = tuple(t[i] + (1 - s[i]) // 2 for i in range(N))
            self.children.append((g, tau))
        # sanity: children partition 2T
        allc = []
        for ch in self.children: allc += self.cubes(ch)
        big = {tuple(2 * x + 1 for x in c) for c in scaled_cells(cells, scale)}
        assert len(allc) == len(set(allc)) == len(big) and set(allc) == big, "children do not partition kT"

    # -- geometry
    def cubes(self, pose):
        g, t = pose; tt = vscale(2, t)
        return [vadd(apply(g, c), tt) for c in self.local_cubes]

    def sub_children(self, G, T, k):
        return [(compose(G, h), vadd(T, vscale(self.scale ** (k - 1), apply(G, tau)))) for (h, tau) in self.children]

    def patch(self, k):
        level = [(self.id, tuple([0] * self.N))]
        for lev in range(k, 0, -1):
            nxt = []
            for (G, T) in level: nxt += self.sub_children(G, T, lev)
            level = nxt
        return level

    def rel(self, A, B):
        gA, tA = A; gB, tB = B; gi = inverse(gA)
        return (compose(gi, gB), apply(gi, vsub(tB, tA)))

    def contacts(self, tiles):
        owner = {}
        for idx, T in enumerate(tiles):
            for c in self.cubes(T): owner[c] = idx
        poses = set()
        for idx, T in enumerate(tiles):
            for c in self.cubes(T):
                for d in self.dirs:
                    j = owner.get(vadd(c, d))
                    if j is not None and j != idx: poses.add(self.rel(T, tiles[j]))
        return poses

    def panel_triples(self, pose):
        g, t = pose; gi = inverse(g); Bc = set(self.cubes(pose)); trip = set()
        for (c, d) in self.panels:
            nb = vadd(c, d)
            if nb in Bc:
                lc = apply(gi, vsub(nb, vscale(2, t))); ld = apply(gi, tuple(-x for x in d))
                trip.add((self.panel_index[(c, d)], self.panel_index[(lc, ld)], g))
        return trip

    def overlaps(self, pose):
        return not set(self.cubes(pose)).isdisjoint(self.local_set)

    def all_touching_poses(self):
        base = self.local_set; out = []
        for g in self.G:
            img = [apply(g, c) for c in self.local_cubes]
            rs = []
            for i in range(self.N):
                lo = min(c[i] for c in img); hi = max(c[i] for c in img)
                tmin = -((hi + 1) // 2) - 1; tmax = (self.hi[i] + 2 - lo) // 2 + 1
                rs.append(range(tmin, tmax + 1))
            for t in itertools.product(*rs):
                tt = vscale(2, t); cs = [vadd(c, tt) for c in img]
                if any(c in base for c in cs): continue
                cset = set(cs)
                if any(vadd(c, d) in cset for c in self.local_cubes for d in self.dirs):
                    out.append((g, t))
        return out

# ---------------------------------------------------------------- dissections
def orientations(cells, G):
    seen = {}
    for g in G:
        o = normalize_cells([apply(g, c) for c in cells])
        seen.setdefault(o, g)
    return seen

def dissections(cells, G=None, limit=10 ** 6, k=2):
    """all exact covers of 2*cells by congruent copies; each cover is a tuple of (g, t) in cell coords,
    t chosen so that g(cells)+t is the placed copy (g = a representative modulo the stabiliser)."""
    cells = normalize_cells(cells); N = len(next(iter(cells)))
    if G is None: G = make_group(N)
    oris = orientations(cells, G)
    tgt = scaled_cells(cells, k); T = set(tgt); bb = [max(c[i] for c in tgt) for i in range(N)]
    pls = []
    for o, g in oris.items():
        img = [apply(g, c) for c in cells]
        m = [min(c[i] for c in img) for i in range(N)]
        mb = [max(c[i] for c in o) for i in range(N)]
        for t in itertools.product(*[range(bb[i] - mb[i] + 1) for i in range(N)]):
            cs = frozenset(tuple(c[i] + t[i] for i in range(N)) for c in o)
            if cs <= T:
                pls.append((cs, g, tuple(t[i] - m[i] for i in range(N))))   # g(cells) + t' = cs
    cell_to_pl = {c: [i for i, (cs, _, _) in enumerate(pls) if c in cs] for c in tgt}
    sols = []
    def rec(rem, chosen):
        if len(sols) >= limit: return
        if not rem: sols.append(tuple(chosen)); return
        c = min(rem, key=lambda c: sum(1 for i in cell_to_pl[c] if pls[i][0] <= rem))
        for i in cell_to_pl[c]:
            if pls[i][0] <= rem:
                chosen.append(i); rec(rem - pls[i][0], chosen); chosen.pop()
    rec(frozenset(tgt), [])
    return [[(pls[i][1], pls[i][2]) for i in sol] for sol in sols]

def dissection_classes(cells, G=None, k=2, limit=10 ** 6):
    """dissections modulo the symmetry group of the scaled tile (acting on the set of child cell-sets)."""
    cells = normalize_cells(cells); N = len(next(iter(cells)))
    if G is None: G = make_group(N)
    tgt = scaled_cells(cells, k)
    syms = stabiliser(tgt, G)          # (g, t) with g(tgt)+t = tgt
    ds = dissections(cells, G, limit=limit, k=k)
    def child_sets(d):
        return frozenset(frozenset(vadd(apply(g, c), t) for c in cells) for (g, t) in d)
    classes = {}
    for d in ds:
        cs = child_sets(d)
        key = min(tuple(sorted(tuple(sorted(frozenset(vadd(apply(g, c), t) for c in S))) for S in cs)) for (g, t) in syms)
        classes.setdefault(key, d)
    return list(classes.values()), len(ds), len(syms)

def frame_variants(cells, dissection, G=None):
    """generator over all frame choices: child i's map g_i may be replaced by g_i*s for s in Stab(T)
    (the copy occupies the same cells). Child 0 is kept fixed (global conjugation)."""
    cells = normalize_cells(cells); N = len(next(iter(cells)))
    if G is None: G = make_group(N)
    stab = stabiliser(cells, G)       # (s, u): s(cells)+u = cells
    k = len(dissection)
    for choice in itertools.product(range(len(stab)), repeat=k - 1):
        d = [dissection[0]]
        for i, ci in enumerate(choice, start=1):
            g, t = dissection[i]; s, u = stab[ci]
            # g(s(c)+u) + t' = g(c') + t  for c' in cells  =>  new map g*s with translation t + g(u)
            d.append((compose(g, s), vadd(t, apply(g, u))))
        yield d

def homochiral(dissection):
    return len({det(g) for g, _ in dissection}) == 1

if __name__ == '__main__':
    # smoke test: the chair reproduces chairN's unique dissection
    N = 3
    chair = frozenset(u for u in itertools.product((0, 1), repeat=N) if u != (1, 1, 1))
    cl, n, ns = dissection_classes(chair)
    print('chair: dissections', n, 'classes', len(cl), 'sym(2T)', ns)
    pt = PolyTile(chair, cl[0], 'chair')
    print('children', pt.children)
    print('touching poses', len(pt.all_touching_poses()), 'panels', len(pt.panels))
