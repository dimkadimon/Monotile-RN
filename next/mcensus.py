import sys, time, pickle
sys.path.insert(0,'/home/user/research'); sys.path.insert(0,'/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from chairN import make_group
from polytile import normalize_cells, dissection_classes, stabiliser
from census3 import SHAPES
from mcert import run_marked
from cert3 import log
name=sys.argv[1]; maxc=int(sys.argv[2]) if len(sys.argv)>2 else 10**9; ML=int(sys.argv[3]) if len(sys.argv)>3 else 4; START=int(sys.argv[4]) if len(sys.argv)>4 else 0
cells=normalize_cells(frozenset(SHAPES[name])); G=make_group(3)
classes,ntot,nsym=dissection_classes(cells,G)
log(f"### {name}: {ntot} dissections -> {len(classes)} classes (frame variants irrelevant under full index marking)")
res=[]
for di,d in list(enumerate(classes))[START:maxc]:
    try:
        r=run_marked(cells,d,f'{name}_c{di}',maxlevels=ML)
    except Exception as e:
        r=dict(status=f'ERR {e!r}'); log(f'{name}_c{di} ERR {e!r}')
    res.append((di,d,r)); pickle.dump(res,open(f'mcensus_{name}_{START}.pkl','wb'))
from collections import Counter
log(f"### {name} summary: {Counter(r['status'].split('(')[0] for _,_,r in res)}")
