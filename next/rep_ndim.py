"""rep_ndim.py -- test whether a polycube in Z^N is rep-2^N (2x copy tiled by 2^N congruent copies, group B_N)."""
import sys, itertools
def group(N):
    return [(p, s) for p in itertools.permutations(range(N)) for s in itertools.product((1,-1), repeat=N)]
def apply(g, c):
    p, s = g; return tuple(s[i]*c[p[i]] for i in range(len(p)))
def normalize(cells):
    N = len(next(iter(cells))); m = [min(c[i] for c in cells) for i in range(N)]
    return frozenset(tuple(c[i]-m[i] for i in range(N)) for c in cells)
def scaled(cells):
    N = len(next(iter(cells)))
    return frozenset(tuple(2*c[i]+d[i] for i in range(N)) for c in cells for d in itertools.product((0,1), repeat=N))
def count_covers(cells, limit=1, want_first=True):
    N = len(next(iter(cells))); G = group(N)
    oris = {}
    for g in G:
        o = normalize([apply(g, c) for c in cells]); oris.setdefault(o, g)
    tgt = scaled(cells); T = set(tgt); bb = [max(c[i] for c in tgt) for i in range(N)]
    pls = []
    for o, g in oris.items():
        mb = [max(c[i] for c in o) for i in range(N)]
        for t in itertools.product(*[range(bb[i]-mb[i]+1) for i in range(N)]):
            cs = frozenset(tuple(c[i]+t[i] for i in range(N)) for c in o)
            if cs <= T: pls.append((cs, g, t))
    cell_to_pl = {c: [i for i,(cs,_,_) in enumerate(pls) if c in cs] for c in tgt}
    sols = []
    def rec(rem, chosen):
        if len(sols) >= limit: return
        if not rem: sols.append(list(chosen)); return
        c = min(rem, key=lambda c: sum(1 for i in cell_to_pl[c] if pls[i][0] <= rem))
        for i in cell_to_pl[c]:
            if pls[i][0] <= rem:
                chosen.append(i); rec(rem - pls[i][0], chosen); chosen.pop()
    rec(frozenset(tgt), [])
    return len(oris), len(pls), sols, pls
if __name__ == '__main__':
    N = int(sys.argv[1]); which = sys.argv[2]
    if which == 'tripod': cells = frozenset([tuple([0]*N)] + [tuple(1 if j==i else 0 for j in range(N)) for i in range(N)])
    elif which == 'chair': cells = frozenset(u for u in itertools.product((0,1), repeat=N) if u != tuple([1]*N))
    elif which == 'screw':  # (0,0,0),(0,0,1),(0,1,0),(1,0,1) generalised: a path 0,e_N, e_N+e_{N-1}... use 3D-literal lift: chain e_i steps
        cells = frozenset([tuple([0]*N)] + [tuple(1 if j>=N-k else 0 for j in range(N)) for k in range(1,N)] + [tuple(1 if (j==0 or j==N-1) else 0 for j in range(N))]) if N>3 else frozenset([(0,0,0),(0,0,1),(0,1,0),(1,0,1)])
    elif which == 'Lprism':  # L-tromino x [0,1]^{N-2}
        cells = frozenset((a,b)+rest for (a,b) in [(0,0),(0,1),(1,0)] for rest in itertools.product((0,),repeat=N-2))
    else: cells = frozenset(eval(which))
    lim = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    no, npl, sols, pls = count_covers(cells, limit=lim)
    print(f"N={N} {which}: cells={len(cells)} orientations={no} placements={npl} dissections found={len(sols)} (limit {lim})")
    if sols:
        gs = [pls[i][1] for i in sols[0]]
        print("  first dissection child orientations:", gs)
