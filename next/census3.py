"""
census3.py -- run the frame-marking certificate over all dissection classes (and frame variants) of a
rep-8 polycube.   usage: python3 census3.py NAME [maxframes] [homochiral_only]
"""
import sys, time, pickle, itertools
sys.path.insert(0, '/home/user/research'); sys.path.insert(0, '/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from chairN import det, make_group
from polytile import PolyTile, dissection_classes, frame_variants, stabiliser, normalize_cells, homochiral
from certpoly import run_tile
from cert3 import log

SHAPES = {
    'tripod':   [(0,0,0),(0,0,1),(0,1,0),(1,0,0)],
    'screw':    [(0,0,0),(0,0,1),(0,1,0),(1,0,1)],
    'hex223':   [(0,0,0),(0,0,1),(0,0,2),(0,1,0),(0,1,1),(1,1,0)],
    'hex222':   [(0,0,0),(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1)],
    'Ltri':     [(0,0,0),(0,0,1),(0,1,0)],
    'Ltetra':   [(0,0,0),(0,0,1),(0,0,2),(0,1,0)],
    'Ppenta':   [(0,0,0),(0,0,1),(0,0,2),(0,1,0),(0,1,1)],
    'Lhexa':    [(0,0,0),(0,0,1),(0,0,2),(0,0,3),(0,1,0),(0,1,1)],
    'chair':    [u for u in itertools.product((0,1), repeat=3) if u != (1,1,1)],
}

def main():
    name = sys.argv[1]; maxframes = int(sys.argv[2]) if len(sys.argv) > 2 else 10 ** 9
    homo_only = len(sys.argv) > 3 and sys.argv[3] == 'homo'
    live = len(sys.argv) > 4 and sys.argv[4] == 'live'
    cells = normalize_cells(frozenset(SHAPES[name]))
    G = make_group(3)
    stab = stabiliser(cells, G)
    classes, ntotal, nsym = dissection_classes(cells, G)
    log(f"### {name}: {len(cells)} cells, |Stab(T)|={len(stab)}, |Sym(2T)|={nsym}, dissections {ntotal} -> {len(classes)} classes")
    results = []; t0 = time.time()
    for di, d in enumerate(classes):
        nvar = len(stab) ** (len(d) - 1)
        chir = 'homochiral' if homochiral(d) else 'mixed'
        if homo_only and not homochiral(d): continue
        log(f"--- dissection class {di}/{len(classes)} ({chir}); frame variants {nvar}")
        import random
        variants = frame_variants(cells, d, G)
        if nvar > maxframes:   # random sample of frame choices (seeded)
            rnd = random.Random(1000 + di)
            k = len(d)
            def sampled():
                from chairN import compose, apply, vadd
                for _ in range(maxframes):
                    fd = [d[0]]
                    for i in range(1, k):
                        g, t = d[i]; s_, u = stab[rnd.randrange(len(stab))]
                        fd.append((compose(g, s_), vadd(t, apply(g, u))))
                    yield fd
            variants = sampled()
        for fi, fd in enumerate(variants):
            if fi >= maxframes: break
            tag = f"{name}_d{di}_f{fi}"
            try:
                r = run_tile(cells, fd, tag, quiet=True, save=False, live_filter=live)
            except AssertionError as e:
                log(f"  {tag}: assertion {e}"); continue
            r['chirality'] = chir; r['class'] = di; r['frame'] = fi
            results.append(r)
            s = r['closure_status']
            if s.startswith('closed'):
                log(f"  {tag}: {s} hist={r['closure_hist']} |P*|={len(r['Pstar'])} tight={r.get('tight')} exact={r.get('coarsen_exact')} "
                    f"claim1={r.get('claim1')} CERTIFIED={r.get('certified')}")
            elif fi < 3 or fi % 50 == 0:
                log(f"  {tag}: {s} hist={r['closure_hist']}")
        pickle.dump(results, open(f'census_{name}{"_live" if live else ""}.pkl', 'wb'))
    # summary
    from collections import Counter
    cnt = Counter(r['closure_status'].split('(')[0] for r in results)
    cert = [r['tag'] for r in results if r.get('certified')]
    log(f"### {name} summary: runs={len(results)} statuses={dict(cnt)} certified={cert} ({time.time()-t0:.0f}s)")

if __name__ == '__main__':
    main()
