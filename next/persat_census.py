import sys, time, pickle
sys.path.insert(0,'/home/user/research'); sys.path.insert(0,'/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from persat import *
from census3 import SHAPES
name=sys.argv[1]; maxc=int(sys.argv[2]); pers=[tuple(int(x) for x in p.split('x')) for p in sys.argv[3:]]
cells=normalize_cells(frozenset(SHAPES[name])); G=make_group(3)
classes,ntot,_=dissection_classes(cells,G)
log(f"### {name}: {ntot} dissections, {len(classes)} classes; k=0; tori {pers}")
out=[]
for ci,D in enumerate(classes[:maxc]):
    H=Hier(cells,D); P=H.adjacency_language(); L=H.language(0); found=None
    for per in pers:
        t0=time.time(); m=torus_search(H,L,0,per,1)
        if m: 
            lv,ok=desubstitute(H,m[0],0,per); found=(per,lv,ok,m[0]); log(f"class {ci}: |P|={len(P)} periodic P-tiling on {per} ({time.time()-t0:.0f}s); desubstitution {'ok' if ok else 'FAILS'} at level {lv}"); break
    if not found: log(f"class {ci}: |P|={len(P)} no periodic P-tiling on tested tori")
    out.append((ci,D,found)); pickle.dump(out,open(f'persat_census_{name}.pkl','wb'))
log(f"### {name}: periodic found for {sum(1 for _,_,f in out if f)}/{len(out)} classes")
