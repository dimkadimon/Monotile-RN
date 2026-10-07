"""rep_census_gen.py N K MAXC CAP [MINC] -- free N-D polycubes <= MAXC cells; test K-dissection (KT by K^N copies);
count dissections (capped) and classes modulo Sym(KT)."""
import sys, itertools, time, pickle
sys.path.insert(0,'/home/user/research/next'); sys.path.insert(0,'/home/user/research')
from polytile import dissections, dissection_classes, normalize_cells
from chairN import make_group, apply
N=int(sys.argv[1]); K=int(sys.argv[2]); MAXC=int(sys.argv[3]); CAP=int(sys.argv[4]); MINC=int(sys.argv[5]) if len(sys.argv)>5 else 2
G=make_group(N)
def canon(cells): return min(tuple(sorted(normalize_cells([apply(g,c) for c in cells]))) for g in G)
shapes={1:{canon([tuple([0]*N)])}}
for k in range(2,MAXC+1):
    new=set()
    for s in shapes[k-1]:
        S=set(s)
        for c in S:
            for d in range(N):
                for e in (-1,1):
                    nc=list(c); nc[d]+=e; nc=tuple(nc)
                    if nc not in S: new.add(canon(S|{nc}))
    shapes[k]=new
print(f"N={N} K={K} free polycubes by size: {[len(shapes[k]) for k in range(1,MAXC+1)]}", flush=True)
res=[]
for k in range(MINC,MAXC+1):
    nrep=0; t0=time.time()
    for s in sorted(shapes[k]):
        cells=frozenset(s)
        ds=dissections(cells,G,limit=CAP,k=K)
        if not ds: continue
        nrep+=1
        bb=tuple(max(c[i] for c in cells)+1 for i in range(N))
        nc=None
        if len(ds)<CAP:
            classes,ntot,nsym=dissection_classes(cells,G,k=K,limit=CAP); nc=len(classes)
        res.append((k,sorted(cells),bb,len(ds),nc))
        print(f"size {k} cells={sorted(cells)} bbox={bb} dissections={len(ds)}{'+' if len(ds)>=CAP else ''} classes={nc}", flush=True)
    print(f"-- size {k}: {nrep} rep-{K**N} shapes of {len(shapes[k])} ({time.time()-t0:.0f}s)", flush=True)
    pickle.dump(res,open(f'rep_census_N{N}_K{K}_{MINC}_{MAXC}.pkl','wb'))
