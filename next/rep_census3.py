"""rep_census3.py -- 3D free polycubes <= MAXC cells: rep-8 test, dissection count (capped) and classes mod Sym(2T).
usage: python3 rep_census3.py MAXC CAP [MINC]"""
import sys, itertools, time, pickle
sys.path.insert(0,'/home/user/research/next'); sys.path.insert(0,'/home/user/research')
from polytile import dissections, dissection_classes, normalize_cells
from chairN import make_group, apply
N=3; MAXC=int(sys.argv[1]); CAP=int(sys.argv[2]); MINC=int(sys.argv[3]) if len(sys.argv)>3 else 2
G=make_group(N)
def canon(cells): return min(tuple(sorted(normalize_cells([apply(g,c) for c in cells]))) for g in G)
shapes={1:{canon([(0,0,0)])}}
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
print(f"free polycubes by size: {[len(shapes[k]) for k in range(1,MAXC+1)]}", flush=True)
res=[]
for k in range(MINC,MAXC+1):
    nrep=0; t0=time.time()
    for s in sorted(shapes[k]):
        cells=frozenset(s)
        ds=dissections(cells,G,limit=CAP)
        if not ds: continue
        nrep+=1
        bb=tuple(max(c[i] for c in cells)+1 for i in range(N))
        if len(ds)<CAP:
            classes,ntot,nsym=dissection_classes(cells,G,limit=CAP); nc=len(classes)
        else: nc=None
        res.append((k,sorted(cells),bb,len(ds),nc))
        print(f"size {k} cells={sorted(cells)} bbox={bb} dissections={len(ds)}{'+' if len(ds)>=CAP else ''} classes={nc}", flush=True)
    print(f"-- size {k}: {nrep} rep-8 shapes of {len(shapes[k])} ({time.time()-t0:.0f}s)", flush=True)
    pickle.dump(res,open(f'rep_census3_{MINC}_{MAXC}.pkl','wb'))
