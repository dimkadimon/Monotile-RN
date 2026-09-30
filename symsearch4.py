"""Enumerate homochiral frame assignments for N=4 invariant under conjugation by a group Gamma of coordinate
permutations (pi_{gamma v} = gamma pi_v gamma^-1), evaluate the closure proxy (number of pairwise-admissible
supertile poses at iteration 0), and store all results."""
import itertools, pickle, sys, time
from chairN import *
from cert3 import language, Rules
from cert6 import supertile_poses, perm_sign

import os
N = int(os.environ.get('NDIM', '4'))
perms = list(itertools.permutations(range(N)))
def pcompose(a, b): return tuple(a[b[i]] for i in range(N))      # (a b)(i) = a[b[i]]  as functions
def pinv(a):
    r = [0] * N
    for i, x in enumerate(a): r[x] = i
    return tuple(r)
def gen_group(gens):
    G = {tuple(range(N))}; frontier = list(G)
    while frontier:
        g = frontier.pop(); 
        for s in gens:
            h = pcompose(s, g)
            if h not in G: G.add(h); frontier.append(h)
    return sorted(G)
def act_v(g, v):  # permute coordinates of v: (g.v)_{g(i)} = v_i
    r = [0] * N
    for i in range(N): r[g[i]] = v[i]
    return tuple(r)
def conj(g, pi):  # frame permutations are used via apply(perm,x)_i = x[perm[i]]  (i.e. matrix of perm^-1);
    # conjugating the linear map: need pi' with apply(pi',x) = g apply(pi, g^-1 x). Work it out numerically.
    def lin(p, x): return tuple(x[p[i]] for i in range(N))
    gi = pinv(g)
    cols = []
    e = [tuple(1 if j == i else 0 for j in range(N)) for i in range(N)]
    # find p' such that lin(p', x) = permute_coords(g, lin(p, permute_coords(gi, x)))
    def pc(g, x):
        r = [0] * N
        for i in range(N): r[g[i]] = x[i]
        return tuple(r)
    for p2 in perms:
        if all(lin(p2, x) == pc(g, lin(pi, pc(gi, x))) for x in e): return p2
    raise RuntimeError

def enumerate_assignments(gens):
    Gam = gen_group(gens)
    U = [u for u in itertools.product((0, 1), repeat=N) if any(c == 0 for c in u)]
    orbits = []; seen = set()
    for v in U:
        if v in seen: continue
        orb = sorted(set(act_v(g, v) for g in Gam)); seen |= set(orb); orbits.append(orb)
    opts = []
    for orb in orbits:
        v = orb[0]
        stab = [g for g in Gam if act_v(g, v) == v]
        par = 1 if sum(v) % 2 == 0 else -1
        ok = [p for p in perms if perm_sign(p) == par and all(conj(g, p) == p for g in stab)]
        if v == (0,) * N: ok = [tuple(range(N))]
        opts.append(ok)
    total = 1
    for o in opts: total *= len(o)
    return Gam, orbits, opts, total

def build(orbits, opts, choice, Gam):
    fr = {'0': tuple(range(N))}
    for orb, p in zip(orbits, choice):
        v0 = orb[0]
        for g in Gam:
            fr[act_v(g, v0)] = conj(g, p)
    return fr

def evaluate(fr):
    ch = Chair(N, fr)
    P, k = language(ch, verbose=False)
    R = Rules(ch, ch.local_cubes, P)
    cand = supertile_poses(ch, R)
    odd_ = sum(1 for q in cand if any(x % 2 for x in q[1]))
    new = set((g, tuple(x // 2 for x in t)) for g, t in cand if not any(x % 2 for x in t)) - set(P)
    return (len(cand), odd_, len(new), len(P))

if __name__ == "__main__":
    gens = eval(sys.argv[1]); tag = sys.argv[2]
    lo = int(sys.argv[3]) if len(sys.argv) > 3 else 0; hi = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9
    Gam, orbits, opts, total = enumerate_assignments(gens)
    print(f"|Gamma|={len(Gam)} orbits={[len(o) for o in orbits]} options={[len(o) for o in opts]} total={total}", flush=True)
    res = []; t0 = time.time(); best = None
    for n, choice in enumerate(itertools.product(*opts)):
        if n < lo or n >= hi: continue
        fr = build(orbits, opts, choice, Gam)
        r = evaluate(fr); res.append((r, fr))
        if best is None or r < best: best = r; print(f"  {n}: new best {r} {fr}", flush=True)
        if n % 50 == 0:
            print(f"  {n}/{total} {time.time()-t0:.0f}s", flush=True)
            pickle.dump(res, open(f'symsearch4_{tag}.pkl', 'wb'))
    pickle.dump(res, open(f'symsearch4_{tag}.pkl', 'wb'))
    print("done; best", best)
