"""rep_census_nd.py -- enumerate free polycubes in Z^N with <= MAXC cells, test rep-2^N, count dissections (capped)."""
import sys, itertools, time
sys.path.insert(0,'/home/user/research/next'); sys.path.insert(0,'/home/user/research')
from polytile import dissections, normalize_cells
from chairN import make_group, apply
N = int(sys.argv[1]); MAXC = int(sys.argv[2]); CAP = int(sys.argv[3]) if len(sys.argv)>3 else 2000
G = make_group(N)
def canon(cells): return min(tuple(sorted(normalize_cells([apply(g,c) for c in cells]))) for g in G)
shapes = {1: {canon([tuple([0]*N)])}}
for k in range(2, MAXC+1):
    new=set()
    for s in shapes[k-1]:
        S=set(s)
        for c in S:
            for d in range(N):
                for e in (-1,1):
                    nc=list(c); nc[d]+=e; nc=tuple(nc)
                    if nc not in S: new.add(canon(S|{nc}))
    shapes[k]=new
print(f"N={N} free polycubes by size: {[len(shapes[k]) for k in range(1,MAXC+1)]}", flush=True)
for k in range(2, MAXC+1):
    for s in sorted(shapes[k]):
        cells=frozenset(s); t0=time.time()
        ds = dissections(cells, G, limit=CAP)
        if not ds: continue
        bb = tuple(max(c[i] for c in cells)+1 for i in range(N))
        ntr = sum(1 for d in ds if all(g==(tuple(range(N)),tuple([1]*N)) for g,_ in d))
        print(f"size {k} cells={sorted(cells)} bbox={bb} prism_dims={sum(1 for x in bb if x==1)} dissections={len(ds)}{'+' if len(ds)>=CAP else ''} translation_only={ntr} ({time.time()-t0:.0f}s)", flush=True)
