"""boundary-consistency objective: for each child panel on the supertile boundary compare its F-partner set
(in absolute frame) with that of the tile panel it covers."""
import itertools
from chairN import *
from cert3 import language, Rules, facet_rules
from cert6 import supertile_poses

def boundary_violations(ch, P, F=None, detail=False):
    N = ch.N
    if F is None: F = facet_rules(ch, P)
    part = {}
    for (a, b, g) in F: part.setdefault(a, set()).add((b, g))
    owner = {}
    for (h, tau) in ch.children:
        for c in ch.local_cubes:
            owner[vadd(apply(h, c), vscale(2, tau))] = (h, tau)
    bad = 0; tot = 0; bad_children = {}
    for (c, d) in ch.panels:
        u = tuple((x - 1) // 2 for x in c)
        i = [k for k in range(N) if d[k] != 0][0]; wi = 1 if d[i] > 0 else 0
        aT = ch.panel_index[(c, d)]
        ST = part.get(aT, set())
        for w in itertools.product((0, 1), repeat=N):
            if w[i] != wi: continue
            x = tuple(4 * u[k] + 2 * w[k] + 1 for k in range(N))
            h, tau = owner[x]
            hi = inverse(h)
            lc = apply(hi, vsub(x, vscale(2, tau))); ld = apply(hi, d)
            aC = ch.panel_index[(lc, ld)]
            SC = set((b, compose(h, g)) for (b, g) in part.get(aC, set()))
            tot += 1
            if SC != ST:
                bad += 1; key = (h, tau); bad_children[key] = bad_children.get(key, 0) + 1
    return (bad, tot, bad_children) if detail else bad

def score(ch):
    P, k = language(ch, verbose=False)
    return boundary_violations(ch, P), len(P)

if __name__ == "__main__":
    import pickle
    res = pickle.load(open('scan3.pkl', 'rb'))
    import collections
    tab = collections.defaultdict(list)
    for fr, p, cand, odd, new in res[:2187]:
        ch = Chair(3, fr)
        P, k = language(ch, verbose=False)
        b = boundary_violations(ch, P)
        tab[cand].append(b)
    for cand in sorted(tab):
        print(cand, collections.Counter(tab[cand]))
