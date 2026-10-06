"""
rep8_polycubes.py -- enumerate small polycubes and test which are rep-8 rep-tiles
(the 2x scaled copy is tiled by 8 congruent copies), with or without reflections.

For every rep-8 polycube we also record every dissection (up to the symmetry of the
scaled tile) and whether the dissection NEEDS rotated children (if all children are
translates the substitution is a lattice substitution and can only produce periodic
hierarchies -- useless for aperiodicity).

usage: python3 rep8_polycubes.py MAXCELLS [chiral]
"""
import sys, itertools
from collections import deque

MAXC = int(sys.argv[1]) if len(sys.argv) > 1 else 6
CHIRAL = len(sys.argv) > 2 and sys.argv[2] == 'chiral'   # only proper rotations

# ---------------------------------------------------------------- symmetry group of Z^3
def group(proper):
    els = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            # parity of permutation
            p = sum(1 for i in range(3) for j in range(i) if perm[j] > perm[i]) % 2
            d = (-1) ** p
            for s in signs: d *= s
            if proper and d != 1: continue
            els.append((perm, signs))
    return els
G24 = group(True); G48 = group(False)

def apply(g, c):
    perm, signs = g
    return tuple(signs[i] * c[perm[i]] for i in range(3))

def normalize(cells):
    m = [min(c[i] for c in cells) for i in range(3)]
    return frozenset(tuple(c[i] - m[i] for i in range(3)) for c in cells)

def canon(cells, G):
    return min((tuple(sorted(normalize([apply(g, c) for c in cells]))) for g in G))

# ---------------------------------------------------------------- polycube enumeration (fixed up to G48)
def polycubes(n):
    shapes = {1: {canon([(0, 0, 0)], G48)}}
    for k in range(2, n + 1):
        new = set()
        for s in shapes[k - 1]:
            S = set(s)
            for c in S:
                for d in range(3):
                    for e in (-1, 1):
                        nc = list(c); nc[d] += e; nc = tuple(nc)
                        if nc not in S:
                            new.add(canon(S | {nc}, G48))
        shapes[k] = new
    return shapes

# ---------------------------------------------------------------- rep-8 test by exact cover
def scaled(cells):
    return frozenset((2 * x + dx, 2 * y + dy, 2 * z + dz) for (x, y, z) in cells
                     for dx in (0, 1) for dy in (0, 1) for dz in (0, 1))

def orientations(cells, G):
    seen = {}
    for g in G:
        o = normalize([apply(g, c) for c in cells])
        if o not in seen: seen[o] = g
    return seen  # oriented shape -> one group element producing it

def placements(target, oris):
    """all placements (oriented shape, translation) inside target; returns list of (frozenset cells, (ori, t))"""
    res = []
    T = set(target)
    bb = [max(c[i] for c in target) for i in range(3)]
    for o, g in oris.items():
        mb = [max(c[i] for c in o) for i in range(3)]
        for tx in range(bb[0] - mb[0] + 1):
            for ty in range(bb[1] - mb[1] + 1):
                for tz in range(bb[2] - mb[2] + 1):
                    cs = frozenset((x + tx, y + ty, z + tz) for (x, y, z) in o)
                    if cs <= T: res.append((cs, (g, (tx, ty, tz))))
    return res

def exact_covers(target, pls, limit=10000):
    """enumerate exact covers (list of placement indices) of target by placements; simple DLX-like backtracking"""
    cell_to_pl = {c: [] for c in target}
    for i, (cs, _) in enumerate(pls):
        for c in cs: cell_to_pl[c].append(i)
    sols = []
    def rec(remaining, chosen):
        if len(sols) >= limit: return
        if not remaining:
            sols.append(list(chosen)); return
        # choose the cell with fewest options
        c = min(remaining, key=lambda c: sum(1 for i in cell_to_pl[c] if pls[i][0] <= remaining))
        for i in cell_to_pl[c]:
            cs = pls[i][0]
            if cs <= remaining:
                chosen.append(i); rec(remaining - cs, chosen); chosen.pop()
    rec(frozenset(target), [])
    return sols

def main():
    G = G24 if CHIRAL else G48
    shapes = polycubes(MAXC)
    print(f"# polycubes by size: {[len(shapes[k]) for k in range(1, MAXC+1)]}  (group {'B3+ (24)' if CHIRAL else 'B3 (48)'})")
    for k in range(1, MAXC + 1):
        for s in sorted(shapes[k]):
            cells = frozenset(s)
            tgt = scaled(cells)
            oris = orientations(cells, G)
            pls = placements(tgt, oris)
            sols = exact_covers(tgt, pls, limit=100000)
            if not sols: continue
            # classify: does some dissection use only translates? does some use non-translates?
            ids = []
            for sol in sols:
                gs = [pls[i][1][0] for i in sol]
                ids.append(all(g == ((0, 1, 2), (1, 1, 1)) for g in gs))
            n_tr = sum(ids); n_rot = len(sols) - n_tr
            # bounding box and "is a prism" (one coordinate constant)
            bb = tuple(max(c[i] for c in cells) + 1 for i in range(3))
            prism = min(bb) == 1
            print(f"size {k:2d}  cells={sorted(cells)}  bbox={bb}{'  (prism of a polyomino)' if prism else ''}")
            print(f"         rep-8 dissections (ordered placements): {len(sols)}   translation-only: {n_tr}   with rotated children: {n_rot}   orientations used: {len(oris)}")
    print("done")

if __name__ == '__main__':
    main()
