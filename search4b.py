import itertools, pickle, time, sys, random
from chairN import *
from cert3 import language, Rules
from cert6 import supertile_poses, perm_sign
N = 4
perms = list(itertools.permutations(range(N)))
even = [p for p in perms if perm_sign(p) == 1]; odd = [p for p in perms if perm_sign(p) == -1]
U = [u for u in itertools.product((0, 1), repeat=N) if any(c == 0 for c in u)]
choices = {v: (even if sum(v) % 2 == 0 else odd) for v in U}
ident = tuple(range(N))
cache = {}
def lang_size(fr):
    key = tuple(fr[v] for v in U)
    if key in cache: return cache[key]
    ch = Chair(N, fr); P, k = language(ch, verbose=False)
    cache[key] = len(P); return len(P)
def poses(fr):
    ch = Chair(N, fr); P, k = language(ch, verbose=False)
    R = Rules(ch, ch.local_cubes, P); cand = supertile_poses(ch, R)
    odd_ = sum(1 for q in cand if any(x % 2 for x in q[1]))
    new = set((g, tuple(x // 2 for x in t)) for g, t in cand if not any(x % 2 for x in t)) - set(P)
    return (len(cand), odd_, len(new), len(P))
def descent(fr):
    best = lang_size(fr); improved = True
    while improved:
        improved = False
        order = list(U); random.shuffle(order)
        for v in order:
            if v == (0,)*N: continue
            cur = fr[v]
            for p in choices[v]:
                if p == cur: continue
                fr[v] = p; r = lang_size(fr)
                if r < best: best = r; cur = p; improved = True
            fr[v] = cur
    return fr, best
if __name__ == "__main__":
    seed = int(sys.argv[1]); random.seed(seed)
    out = open(f'search4b_{seed}.log', 'a'); results = []; seenmin = {}
    for restart in range(100000):
        fr = {'0': ident, (0,)*N: ident}
        for v in U:
            if v != (0,)*N: fr[v] = random.choice(choices[v])
        t0 = time.time(); fr, best = descent(fr)
        key = tuple(fr[v] for v in U)
        if key in seenmin: print(f"restart {restart}: min |P|={best} (seen)", file=out, flush=True); continue
        r = poses(fr); seenmin[key] = r
        print(f"restart {restart}: min |P|={best} poses={r} {time.time()-t0:.0f}s {fr}", file=out, flush=True)
        results.append((r, dict(fr)))
        pickle.dump(results, open(f'search4b_{seed}.pkl', 'wb'))
