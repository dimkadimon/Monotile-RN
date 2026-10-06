import sys, pickle, time
sys.path.insert(0,'/home/user/research'); sys.path.insert(0,'/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from persat import *
from census3 import SHAPES
from mcert import run_marked
cells=normalize_cells(frozenset(SHAPES['plate25'])); G=make_group(3)
classes,_,_=dissection_classes(cells,G)
tori=[(2,4,10),(2,8,10),(4,8,10),(2,10,10),(4,10,10),(2,4,20),(5,4,8),(5,8,8),(4,6,10),(6,6,10)]
for ci in [1,2,5,6,7,9]:
    D=classes[ci]; H=Hier(cells,D); P=H.adjacency_language(); L=H.language(0)
    found=None
    for per in tori:
        t0=time.time(); m=torus_search(H,L,0,per,1)
        if m: found=(per,desub_Z(H,m[0],0,per,maxlevel=3)); log(f"class {ci}: periodic P-tiling on {per}: desub {found[1]} ({time.time()-t0:.0f}s)"); break
    if not found: log(f"class {ci}: no periodic P-tiling on {tori}")
    r=run_marked(cells,D,f'plate25_c{ci}',maxlevels=4)
    log(f"class {ci}: mcert k=1 status {r['status']} certified={r.get('certified')}")
