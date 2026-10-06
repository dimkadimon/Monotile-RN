"""
certpoly.py -- frame-marking certificate (closure, tightness, two-shell enclosure analysis) for an
arbitrary polycube rep-2^N tile with a chosen dissection + frames.  Reuses cert6/cert3 machinery.

run_tile(cells, dissection, tag) -> dict with closure status, |P*|, |F*|, tightness, claim-1 stats, certified flag.
"""
import sys, time, pickle, itertools
sys.path.insert(0, '/home/user/research'); sys.path.insert(0, '/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from chairN import apply, compose, inverse, vadd, vscale, identity
from cert3 import Rules, language, facet_rules, log
from cert6 import HCert, closure, supertile_poses, tightness_H
from polytile import PolyTile

def claim1_generic(C, progress=200, quiet=False):
    """Two-shell enclosure analysis.  For every complete enclosure E of the reference tile T:
       dead (not extendable)  |  unique compatible supertile S fully present (type A)
       |  unique S, some kids missing: there is a kid U in the patch such that every enclosure of U
          either completes S or is dead (type B, 'presence')."""
    T = C.origin; t0 = time.time()
    encl = C.enclosures_given(T, [])
    st = dict(enclosures=len(encl), live=0, dead=0, U_fail=0, Pr_fail=0, typeA=0, typeB=0, Pr_checks=0)
    fails = []
    if not quiet: log(f"  enclosures: {len(encl)}")
    for n, E in enumerate(encl):
        patch = [T] + list(E)
        if not C.live(patch): st['dead'] += 1; continue
        st['live'] += 1
        cs = C.compat(T, patch)
        if len(cs) != 1:
            st['U_fail'] += 1; fails.append(('U', E, [c for c, _, _ in cs])); continue
        c, S, K_all = cs[0]
        missing = [K for K in K_all if K not in patch]
        if not missing: st['typeA'] += 1; continue
        st['typeB'] += 1
        # try each present kid as the pivot U
        ok_any = False
        for U in [K for K in K_all if K in patch and K != T]:
            st['Pr_checks'] += 1
            good = True
            for EU in C.enclosures_given(U, patch):
                joint = list(dict.fromkeys(patch + list(EU)))
                if all(K in joint for K in missing): continue
                if C.live(joint): good = False; break
            if good: ok_any = True; break
        if not ok_any: st['Pr_fail'] += 1; fails.append(('Pr', E))
        if not quiet and n % progress == 0: log(f"    {n}/{len(encl)} {st} ({time.time()-t0:.0f}s)")
    st['limit'] = C.stats['limit']; st['time'] = time.time() - t0
    if not quiet: log("  Claim1:", st)
    return st, fails

def closure_live(ch, H, maxit=12, node_limit=20000, quiet=True):
    """Liveness-filtered coarsening closure.  A pairwise-admissible supertile pose is kept only if the
    two-supertile patch is extendable (tile-level DFS liveness under the current rule).  Dropped poses are
    re-checked against the final rule; the result is a fixpoint P* such that every *live* pairwise-admissible
    supertile pose is aligned and lies in 2P*.  Sound because poses occurring in a tiling are always live."""
    from cert3 import language
    P0, k = language(ch, verbose=False)
    P = set(P0); hist = [len(P0)]
    S0 = [(h, tau) for (h, tau) in ch.children]
    def kids(S):
        G, T = S
        return [(compose(G, h), vadd(T, apply(G, tau))) for (h, tau) in ch.children]
    dropped = set(); n_dead_total = 0
    for it in range(maxit):
        R = Rules(ch, ch.local_cubes, P)
        cand = supertile_poses(ch, R)
        C = HCert(ch, H, P, node_limit=node_limit, use_prop=True)
        live = set(); dead = set()
        for S in cand:
            if S in dropped:              # re-test dropped poses under the current (larger) rule
                pass
            patch = S0 + kids(S)
            if C.live(patch): live.add(S)
            else: dead.add(S)
        if C.stats['limit']: return None, hist, f'limit({C.stats["limit"]})', None
        dropped |= dead; n_dead_total = len(dropped)
        odd = [q for q in live if any(x % 2 for x in q[1])]
        if odd: return None, hist, f'offset({len(odd)})/live{len(live)}/dead{len(dead)}', None
        new = set((g, tuple(x // 2 for x in t)) for g, t in live)
        if not quiet: log(f"  closure_live it{it}: |P|={len(P)} cand={len(cand)} live={len(live)} dead={len(dead)} new={len(new - P)}")
        if new <= P:
            hist.append(len(P))
            return frozenset(P), hist, ('closed' if new == P else 'closed-strict') + f'/dead{len(dead)}', live
        P |= new; hist.append(len(P))
    return None, hist, 'diverged', None

def closure_k(ch, H, maxit=10, quiet=True):
    """coarsening closure for inflation factor k = ch.scale (cert6.closure generalised)."""
    from cert3 import language
    k = ch.scale
    P0, lev = language(ch, verbose=False)
    P = set(P0); hist = [len(P0), len(P)]
    for it in range(maxit):
        R = Rules(ch, ch.local_cubes, P)
        cand = supertile_poses(ch, R)
        odd = [q for q in cand if any(x % k for x in q[1])]
        if odd: return None, hist, f'offset({len(odd)})', cand
        new = set((g, tuple(x // k for x in t)) for g, t in cand)
        if not quiet: log(f"  closure it{it}: |P|={len(P)} supertile poses={len(cand)} new={len(new - P)}")
        if new <= P:
            return frozenset(P), hist, ('closed' if new == P else 'closed-strict'), cand
        P |= new; hist.append(len(P))
    return None, hist, 'diverged', None

def run_tile(cells, dissection, tag, node_limit=50000, quiet=False, do_claim=True, save=True, live_filter=False, scale=2):
    ch = PolyTile(cells, dissection, tag, scale=scale); H = [identity(ch.N)]; t0 = time.time()
    if not quiet: log(f"=== {tag}: N={ch.N} cells={len(ch.cells)} children={len(ch.children)} panels={len(ch.panels)}")
    if live_filter: Pstar, hist, status, cand = closure_live(ch, H, quiet=quiet)
    else: Pstar, hist, status, cand = closure_k(ch, H, quiet=quiet)
    out = dict(tag=tag, cells=sorted(ch.cells), dissection=dissection, closure_status=status, closure_hist=hist)
    if not quiet: log(f"  closure: {status} sizes={hist}")
    if Pstar is None:
        if save: pickle.dump(out, open(f'cp_{tag}.pkl', 'wb'))
        return out
    scaled = set((g, vscale(ch.scale, t)) for g, t in Pstar)
    F = facet_rules(ch, Pstar); ntouch, PF = tightness_H(ch, Pstar, F)
    out.update(Pstar=sorted(Pstar), F=sorted(F), touching=ntouch, tight=(PF == Pstar), coarsen_exact=(cand == scaled))
    if not quiet: log(f"  |P*|={len(Pstar)} |F*|={len(F)} touching={ntouch} |P_F*|={len(PF)} tight={PF == Pstar} exact={cand == scaled} ({time.time()-t0:.0f}s)")
    if PF != Pstar or not do_claim:
        out['PF_extra'] = sorted(PF - Pstar)
        if save: pickle.dump(out, open(f'cp_{tag}.pkl', 'wb'))
        return out
    C = HCert(ch, H, Pstar, node_limit=node_limit, use_prop=True)
    st, fails = claim1_generic(C, quiet=quiet)
    ok = st['U_fail'] == 0 and st['Pr_fail'] == 0 and st['limit'] == 0 and cand == scaled
    out.update(claim1=st, fails=fails[:50], certified=ok)
    if save: pickle.dump(out, open(f'cp_{tag}.pkl', 'wb'))
    if not quiet: log(f"=== {tag}: CERTIFIED={ok}  ({time.time()-t0:.0f}s)")
    return out

if __name__ == '__main__':
    # validation on the chair with the Paper-1 frames: must certify with |P*| = 44
    from chairN import Chair
    fr30 = {'0': (0, 1, 2), (0, 0, 0): (0, 1, 2), (0, 0, 1): (0, 2, 1), (0, 1, 0): (2, 1, 0), (0, 1, 1): (1, 2, 0),
            (1, 0, 0): (1, 0, 2), (1, 0, 1): (0, 1, 2), (1, 1, 0): (2, 0, 1)}
    c = Chair(3, fr30)
    ones = (1, 1, 1)
    dis = []
    for (h, tau) in c.children:
        s = apply(h, ones); dis.append((h, tuple(tau[i] - (1 - s[i]) // 2 for i in range(3))))
    cells = frozenset(u for u in itertools.product((0, 1), repeat=3) if u != (1, 1, 1))
    r = run_tile(cells, dis, 'chair_check', save=False)
    print('certified:', r.get('certified'), '|P*|=', len(r.get('Pstar', [])))
