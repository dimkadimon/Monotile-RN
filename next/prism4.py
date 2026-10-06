import sys, time, itertools, random
sys.path.insert(0,'/home/user/research'); sys.path.insert(0,'/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from chairN import make_group, compose, apply, vadd
from polytile import dissections, stabiliser, normalize_cells, homochiral
from certpoly import run_tile
from cert3 import log
cells = normalize_cells(frozenset([(0,0,0,0),(0,1,0,0),(1,0,0,0)]))
G = make_group(4); stab = stabiliser(cells, G)
t0=time.time(); ds = dissections(cells, G, limit=int(sys.argv[1])); log(f"L-prism 4D: {len(ds)} dissections found ({time.time()-t0:.0f}s), |Stab|={len(stab)}")
nsample = int(sys.argv[2]); rnd = random.Random(7)
from collections import Counter; cnt = Counter()
for di, d in enumerate(ds[:int(sys.argv[3])]):
    for fi in range(nsample):
        fd = [d[0]]
        for i in range(1, len(d)):
            g, t = d[i]; s_, u = stab[rnd.randrange(len(stab))] if fi else stab[0]
            fd.append((compose(g, s_), vadd(t, apply(g, u))))
        r = run_tile(cells, fd, f"Lp4_d{di}_f{fi}", quiet=True, save=False)
        cnt[r['closure_status'].split('(')[0]] += 1
        if r['closure_status'].startswith('closed'): log(f"  d{di} f{fi}: CLOSED hist={r['closure_hist']} |P*|={len(r['Pstar'])} tight={r.get('tight')} cert={r.get('certified')} claim1={r.get('claim1')}")
        elif fi == 0: log(f"  d{di} f0: {r['closure_status']} hist={r['closure_hist']} ({time.time()-t0:.0f}s)")
log("summary", dict(cnt))
