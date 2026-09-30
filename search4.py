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
def evaluate(fr):
    key = tuple(fr[v] for v in U)
    if key in cache: return cache[key]
    ch = Chair(N, fr)
    P, k = language(ch, verbose=False)
    R = Rules(ch, ch.local_cubes, P)
    cand = supertile_poses(ch, R)
    odd_ = sum(1 for q in cand if any(x % 2 for x in q[1]))
    new = set((g, tuple(x // 2 for x in t)) for g, t in cand if not any(x % 2 for x in t)) - set(P)
    r = (len(cand), odd_, len(new), len(P)); cache[key] = r; return r
def descent(fr, log):
    best = evaluate(fr); improved = True
    while improved:
        improved = False
        order = list(U); random.shuffle(order)
        for v in order:
            if v == (0,)*N: continue
            cur = fr[v]
            for p in choices[v]:
                if p == cur: continue
                fr[v] = p; r = evaluate(fr)
                if r < best: best = r; cur = p; improved = True; log(f"  improve {best} v={v} p={p}")
            fr[v] = cur
    return fr, best
if __name__ == "__main__":
    seed = int(sys.argv[1]); random.seed(seed)
    out = open(f'search4_{seed}.log', 'a')
    def log(*a):
        print(*a, file=out, flush=True)
    results = []
    for restart in range(1000):
        fr = {'0': ident, (0,)*N: ident}
        for v in U:
            if v != (0,)*N: fr[v] = random.choice(choices[v])
        t0 = time.time(); log(f"restart {restart} start {evaluate(fr)}")
        fr, best = descent(fr, log)
        log(f"restart {restart} DONE {best} {time.time()-t0:.0f}s {fr}")
        results.append((best, dict(fr)))
        pickle.dump(results, open(f'search4_{seed}.pkl', 'wb'))
