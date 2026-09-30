"""Full coarsening-closure test for all 2187 homochiral N=3 frame assignments with T0 = id (from scan3.pkl)."""
import pickle, itertools, time, collections
from chairN import *
from cert3 import language, Rules
from cert6 import supertile_poses
res = pickle.load(open('scan3.pkl', 'rb'))
out = []; t0 = time.time()

def closure_capped(ch, maxit=4, maxP=300):
    P, k = language(ch, verbose=False); P = set(P); hist = [len(P)]
    for it in range(maxit):
        if len(P) > maxP: return 'diverged(>%d)' % maxP, hist
        R = Rules(ch, ch.local_cubes, P); cand = supertile_poses(ch, R)
        if any(any(x % 2 for x in q[1]) for q in cand): return 'offset', hist + [len(cand)]
        new = set((g, tuple(x // 2 for x in t)) for g, t in cand)
        if new <= P: return ('closed' if new == P else 'closed-strict'), hist
        P |= new; hist.append(len(P))
    return 'diverged', hist

for n, (fr, p, cand, odd, new) in enumerate(res):
    ch = Chair(3, fr)
    status, hist = closure_capped(ch)
    out.append((fr, status, hist))
    if n % 100 == 0:
        print(n, status, hist, f"{time.time()-t0:.0f}s", flush=True); pickle.dump(out, open('closure_scan3.pkl', 'wb'))
pickle.dump(out, open('closure_scan3.pkl', 'wb'))
print(collections.Counter(s for _, s, _ in out))
