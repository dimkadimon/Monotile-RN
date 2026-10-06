"""per dissection class: (1) is the D-hierarchy periodic?  (2) periodic P-tiling (L_0) via SAT on small tori;
(3) is it non-hierarchical (de-substitution in Z^3 fails)?   usage: python3 rigidity_census.py NAME MAXC"""
import sys, time, pickle, itertools
sys.path.insert(0,'/home/user/research'); sys.path.insert(0,'/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from persat import *
from census3 import SHAPES

def hier_period(H, m=5, rng=8, mmax=6):
    S0=(H.id,(0,)*H.N); tiles=H.expand(S0,m,0,lambda l:H.D)
    tset={(g,t) for g,t,_ in tiles}; allc=set()
    for g,t,_ in tiles: allc|=H.cubes(g,t)
    for v in sorted(itertools.product(range(-rng,rng+1),repeat=3), key=lambda v: sum(abs(x) for x in v)):
        if v==(0,0,0): continue
        v2=vscale(2,v); inside=0; match=0; bad=False
        for g,t in tset:
            cb=H.cubes(g,t)
            if all(vadd(c,v2) in allc for c in cb):
                inside+=1
                if (g,vadd(t,v)) in tset: match+=1
                else: bad=True; break
        if not bad and inside>200:
            # confirm at level mmax
            tiles6=H.expand(S0,mmax,0,lambda l:H.D); ts6={(g,t) for g,t,_ in tiles6}; allc6=set()
            for g,t,_ in tiles6: allc6|=H.cubes(g,t)
            ok=all((g,vadd(t,v)) in ts6 for g,t in ts6 if all(vadd(c,v2) in allc6 for c in H.cubes(g,t)))
            if ok: return v
    return None

name=sys.argv[1]; maxc=int(sys.argv[2])
cells=normalize_cells(frozenset(SHAPES[name])); G=make_group(3)
classes,ntot,_=dissection_classes(cells,G)
vol=len(cells)
tori=[p for p in [(4,4,4),(4,4,6),(4,4,8),(4,6,6),(4,4,12),(6,6,6),(6,6,8),(4,8,8),(8,8,8),(4,4,7),(4,7,7),(4,4,14),(7,7,7),(4,7,14),(6,6,9),(4,6,9),(6,6,3),(4,4,5),(4,4,10),(4,5,5),(4,4,15)] if (p[0]*p[1]*p[2])%vol==0]
log(f"### {name}: {ntot} dissections, {len(classes)} classes; tori {tori}")
rows=[]
for ci,D in enumerate(classes[:maxc]):
    t0=time.time(); H=Hier(cells,D); P=H.adjacency_language(); L=H.language(0)
    hv=hier_period(H)
    found=None
    for per in tori:
        m=torus_search(H,L,0,per,1)
        if m:
            lv,ok=desub_Z(H,m[0],0,per,maxlevel=3); found=(per,lv,ok,m[0]); break
    rows.append(dict(ci=ci,D=D,P=len(P),hier_period=hv,periodic=found and found[0],desub=found and (found[1],found[2]),tiles=found and found[3]))
    log(f"class {ci}: |P|={len(P)}; hierarchy periodic: {hv}; periodic P-tiling: {found and found[0]}; de-substitution: {found and ('ok to level %d'%found[1] if found[2] else 'FAILS at level %d'%found[1])} ({time.time()-t0:.0f}s)")
    pickle.dump(rows,open(f'rigidity_census_{name}.pkl','wb'))
n=len(rows); hp=sum(1 for r in rows if r['hier_period']); pf=sum(1 for r in rows if r['periodic']); nh=sum(1 for r in rows if r['periodic'] and not r['desub'][1])
log(f"### {name}: classes {n}; periodic hierarchy {hp}; periodic P-tiling found {pf}; of which provably non-hierarchical {nh}; neither {n-pf}")
