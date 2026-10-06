"""
localsearch.py -- coordinate-descent search over frame assignments of a dissection, minimising the
coarsening growth  score = |new \\ P| summed over the first closure iterations (0 = closed).
usage: python3 localsearch.py NAME restarts sweeps [classes]
"""
import sys, time, random, pickle
sys.path.insert(0, '/home/user/research'); sys.path.insert(0, '/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from chairN import make_group, compose, apply, vadd
from cert3 import Rules, language, log
from cert6 import supertile_poses
from polytile import PolyTile, dissection_classes, stabiliser, normalize_cells, homochiral
from census3 import SHAPES
from certpoly import run_tile

def score(cells, fd, maxit=2):
    """growth score: sum over iterations of |new - P|, +inf-ish penalty for offsets; 0 means closed."""
    try: ch = PolyTile(cells, fd)
    except AssertionError: return 10 ** 9, 'bad'
    P0, _ = language(ch, verbose=False); P = set(P0); tot = 0
    for it in range(maxit):
        R = Rules(ch, ch.local_cubes, P); cand = supertile_poses(ch, R)
        odd = sum(1 for q in cand if any(x % 2 for x in q[1]))
        if odd: return tot + 1000 * odd + len(cand), f'offset{odd}@{it}'
        new = set((g, tuple(x // 2 for x in t)) for g, t in cand) - P
        if not new: return tot, 'closed'
        tot += len(new); P |= new
    return tot, 'open'

def main():
    name = sys.argv[1]; restarts = int(sys.argv[2]); sweeps = int(sys.argv[3])
    cells = normalize_cells(frozenset(SHAPES[name])); G = make_group(3)
    stab = stabiliser(cells, G); classes, n, ns = dissection_classes(cells, G)
    sel = range(len(classes)) if len(sys.argv) < 5 else [int(x) for x in sys.argv[4].split(',')]
    log(f"### localsearch {name}: {len(classes)} classes, |Stab|={len(stab)}")
    best_overall = []
    for di in sel:
        d = classes[di]; k = len(d); rnd = random.Random(di)
        def build(choice):
            fd = [d[0]]
            for i in range(1, k):
                g, t = d[i]; s_, u = stab[choice[i - 1]]
                fd.append((compose(g, s_), vadd(t, apply(g, u))))
            return fd
        best = (10 ** 9, None, None)
        for r in range(restarts):
            choice = [rnd.randrange(len(stab)) for _ in range(k - 1)]
            cur, st = score(cells, build(choice))
            for sw in range(sweeps):
                improved = False
                for i in range(k - 1):
                    for v in range(len(stab)):
                        if v == choice[i]: continue
                        c2 = list(choice); c2[i] = v
                        s2, st2 = score(cells, build(c2))
                        if s2 < cur: cur, st, choice, improved = s2, st2, c2, True
                if cur == 0 or not improved: break
            if cur < best[0]: best = (cur, st, list(choice))
            if cur == 0: break
        log(f"  class {di}: best score {best[0]} ({best[1]}) choice={best[2]}")
        best_overall.append((di, best))
        if best[0] == 0:
            res = run_tile(cells, build(best[2]), f"{name}_ls_d{di}", quiet=False, save=True)
            log(f"  >>> class {di}: CERTIFIED={res.get('certified')} |P*|={len(res.get('Pstar', []))} tight={res.get('tight')} claim1={res.get('claim1')}")
    pickle.dump(best_overall, open(f'localsearch_{name}.pkl', 'wb'))
    log("### done; best scores:", sorted(b[0] for _, b in best_overall)[:10])

if __name__ == '__main__':
    main()
