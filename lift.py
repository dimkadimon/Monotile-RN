import itertools, sys, time
from chairN import *
from cert3 import language, Rules
from cert6 import supertile_poses, perm_sign
fr3 = {(0,0,0):(0,1,2),(0,0,1):(0,2,1),(0,1,0):(2,1,0),(0,1,1):(1,2,0),(1,0,0):(1,0,2),(1,0,1):(0,1,2),(1,1,0):(2,0,1)}
N = 4
def pc(a, b): return tuple(a[b[i]] for i in range(N))
def ext(p): return tuple(p) + (3,)
def tr(a):
    t = list(range(4)); t[a], t[3] = 3, a; return tuple(t)
def evaluate(fr):
    ch = Chair(N, fr); P, k = language(ch, verbose=False)
    R = Rules(ch, ch.local_cubes, P); cand = supertile_poses(ch, R)
    odd_ = sum(1 for q in cand if any(x % 2 for x in q[1]))
    new = set((g, tuple(x // 2 for x in t)) for g, t in cand if not any(x % 2 for x in t)) - set(P)
    return (len(cand), odd_, len(new), len(P))
U = [u for u in itertools.product((0, 1), repeat=N) if any(c == 0 for c in u)]
if __name__ == "__main__":
    for a in range(3):
        for order in ('L', 'R'):
            fr = {'0': (0,1,2,3)}
            for v in U:
                vp = v[:3]
                base = ext(fr3[vp]) if vp != (1,1,1) else (0,1,2,3)
                if v[3] == 0: fr[v] = base
                else: fr[v] = pc(tr(a), base) if order == 'L' else pc(base, tr(a))
            t0 = time.time(); print(a, order, evaluate(fr), f"{time.time()-t0:.0f}s", flush=True)
