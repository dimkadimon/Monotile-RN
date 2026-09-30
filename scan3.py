import itertools, pickle, time, sys
from chairN import *
from cert3 import language, Rules
from cert6 import supertile_poses, perm_sign, symmetrize
N = int(sys.argv[1]); out = sys.argv[2]
perms = list(itertools.permutations(range(N)))
even = [p for p in perms if perm_sign(p) == 1]; odd = [p for p in perms if perm_sign(p) == -1]
U = [u for u in itertools.product((0, 1), repeat=N) if any(c == 0 for c in u)]
choices = [even if sum(v) % 2 == 0 else odd for v in U]   # homochiral: det h_v = +1
res = []; t0 = time.time()
for n, pis in enumerate(itertools.product(*choices)):
    fr = {'0': tuple(range(N))}; fr.update(dict(zip(U, pis)))
    ch = Chair(N, fr)
    P, k = language(ch, verbose=False)
    R = Rules(ch, ch.local_cubes, P)
    cand = supertile_poses(ch, R)
    odd_ = sum(1 for q in cand if any(x % 2 for x in q[1]))
    new = set((g, tuple(x // 2 for x in t)) for g, t in cand if not any(x % 2 for x in t)) - set(P)
    res.append((fr, len(P), len(cand), odd_, len(new)))
    if n % 200 == 0:
        print(n, len(P), len(cand), odd_, len(new), f"{time.time()-t0:.0f}s", flush=True)
        pickle.dump(res, open(out, 'wb'))
pickle.dump(res, open(out, 'wb'))
print("done", len(res))
