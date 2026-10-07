import sys, time
sys.path.insert(0,'/home/user/research'); sys.path.insert(0,'/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from persat import *
from mcert import run_marked
cells=normalize_cells(frozenset([(0,0,0,0),(0,0,0,1),(0,0,0,2),(0,0,1,0),(0,1,0,2),(1,0,1,0)])); G=make_group(4)
classes,ntot,_=dissection_classes(cells,G); log(f"tw6: {ntot} dissections, {len(classes)} classes")
for ci,D in enumerate(classes):
    H=Hier(cells,D); P=H.adjacency_language(); L=H.language(0); log(f"class {ci}: |P|={len(P)} |L0|={len(L)}")
    r=run_marked(cells,D,f'tw6_c{ci}',maxlevels=3); log(f"class {ci}: mcert {r['status']} {r.get('certified')}")
    for per in [(4,4,4,6),(4,4,6,6)]:
        t0=time.time(); m=torus_search(H,L,0,per,1)
        if m: log(f"class {ci}: periodic P-tiling on {per}, desub {desub_Z(H,m[0],0,per,maxlevel=3)} ({time.time()-t0:.0f}s)"); break
        log(f"class {ci}: none on {per} ({time.time()-t0:.0f}s)")
