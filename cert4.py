"""
cert4.py -- two-shell finite certificate (Flicker/Goodman-Strauss style), tile level only.

Setting: N-chair with frame assignment; P = stable contact language (monotone fixed point); rules = "adjacent tiles
must be in a relative pose in P" (equivalently the facet rules F, by tightness).

For a tile T and a complete enclosure E_T (all panels of T covered, pairwise consistent) define
   compat(E_T) = { level-1 supertile placements S containing T as a child such that every child of S is either a
                   tile of E_T ∪ {T} or has cubes disjoint from all tiles of E_T ∪ {T} }.
Certificate conditions, checked over all enclosures E_T that survive a 1-step lookahead and a DFS extension test
("live" enclosures):
  (U)  |compat(E_T)| = 1.                                          -> S(E_T)
  (Pr) Presence: for every child K of S(E_T) that is not in E_T ∪ {T} (only possible when T is an outer child; then all
       children touch the central child U ∈ E_T): for every enclosure E_U of U consistent with E_T ∪ {T} whose joint
       patch is live, K ∈ E_U.   (Then in any tiling all children of S(E_T) are present.)
  (Co) Coarsening: for every neighbour T' ∈ E_T with T' ∉ S(E_T): for every enclosure E_T' of T' consistent with
       E_T ∪ {T} whose joint patch is live, |compat(E_T')| = 1 and pose(S(E_T') rel S(E_T)) ∈ 2P.
Consequences: in any rule-abiding tiling every tile lies in exactly one complete supertile; adjacent supertiles are
in poses 2P, so supertiles (scaled by 1/2) again form a rule-abiding tiling; induction gives the unique hierarchy.
"""
import itertools, time, sys, pickle
from chairN import *
from cert3 import Rules, language, facet_rules, tightness, log

class Cert:
    def __init__(self, ch, P, node_limit=20000):
        self.ch = ch; self.N = ch.N; self.P = P
        self.R = Rules(ch, ch.local_cubes, P)
        self.node_limit = node_limit
        self.origin = (ch.id, tuple([0] * ch.N))
        ch_ = ch.children
        self.kids_local = ch_
        # for child index c: supertile placement (G,T) such that child c sits at the origin pose
        self.sup_for_child = []
        for (h, tau) in ch_:
            G = inverse(h); self.sup_for_child.append((G, tuple(-x for x in apply(G, tau))))
        self.scaledP = frozenset((g, vscale(2, t)) for (g, t) in P)
        self.stats = dict(enclosures=0, live=0, dead=0, limit=0, U_fail=0, Pr_checks=0, Pr_fail=0,
                          Co_checks=0, Co_fail=0, joint_live=0, joint_dead=0)
        self.failures = []

    # ---- geometry helpers
    def kids(self, S):
        G, T = S
        return [(compose(G, h), vadd(T, apply(G, tau))) for (h, tau) in self.kids_local]

    def sup_containing(self, A, c):
        """supertile placement in which tile A is child c."""
        G0, T0 = self.sup_for_child[c]
        gA, tA = A
        return (compose(gA, G0), vadd(tA, apply(gA, T0)))

    def compat(self, A, patch):
        """supertile placements containing tile A as some child, compatible with the patch (set of tiles)."""
        cubes = set()
        for t in patch: cubes |= self.R.cubes(t)
        out = []
        for c in range(len(self.kids_local)):
            S = self.sup_containing(A, c)
            ok = True
            for K in self.kids(S):
                if K in patch: continue
                if not self.R.cubes(K).isdisjoint(cubes): ok = False; break
            if ok: out.append((c, S))
        return out

    def consistent(self, tiles):
        for i in range(len(tiles)):
            for j in range(i + 1, len(tiles)):
                if not self.R.pair_ok(tiles[i], tiles[j]): return False
        return True

    # ---- enclosures of `centre` given fixed tiles
    def enclosures_given(self, centre, fixed):
        R = self.R
        fixed = [t for t in fixed if t != centre]
        cand = [R.absolute(centre, rp) for rp in R.Plist]
        cand = [B for B in cand if B in fixed or all(R.pair_ok(B, f) for f in fixed)]
        M = len(cand)
        comp = [[True] * M for _ in range(M)]
        for i in range(M):
            for j in range(i + 1, M):
                comp[i][j] = comp[j][i] = R.pair_ok(cand[i], cand[j])
        cc = R.cubes(centre)
        panels = [(c, d) for c in cc for d in R.dirs if vadd(c, d) not in cc]
        covers = [frozenset(k for k, (c, d) in enumerate(panels) if vadd(c, d) in R.cubes(cand[i])) for i in range(M)]
        by_panel = [[i for i in range(M) if k in covers[i]] for k in range(len(panels))]
        out = []
        def bt(chosen, covered):
            k = 0
            while k < len(panels) and k in covered: k += 1
            if k == len(panels):
                out.append(tuple(cand[i] for i in chosen)); return
            for i in by_panel[k]:
                if all(comp[i][j] for j in chosen):
                    bt(chosen + [i], covered | covers[i])
        bt([], frozenset())
        return out

    def live(self, patch):
        """patch (list of tiles) is live if it passes lookahead and admits an extension enclosing all its tiles."""
        R = self.R
        occ = {}
        for T in patch:
            for c in R.cubes(T): occ[c] = T
        if not R.lookahead_ok(patch, occ, patch): return False
        ext = R.extend_exists(list(patch), patch, self.node_limit)
        if R.limit_hit:
            self.stats['limit'] += 1
            return True   # undetermined -> treat as live (conservative)
        return ext

    # ---- main
    def run(self, progress=50):
        t0 = time.time()
        T = self.origin
        encl = self.enclosures_given(T, [])
        self.stats['enclosures'] = len(encl)
        log(f"  enclosures of a tile: {len(encl)}")
        for n, E in enumerate(encl):
            patch = [T] + list(E)
            if not self.live(patch):
                self.stats['dead'] += 1; continue
            self.stats['live'] += 1
            cs = self.compat(T, set(patch))
            if len(cs) != 1:
                self.stats['U_fail'] += 1; self.failures.append(('U', E, cs)); continue
            c, S = cs[0]
            K_all = self.kids(S)
            Sset = set(K_all)
            missing = [K for K in K_all if K not in patch]
            # (Pr) presence
            if missing:
                U = K_all[0]  # central child
                assert U in patch, "central child must be in the enclosure of an outer child"
                self.stats['Pr_checks'] += 1
                for EU in self.enclosures_given(U, patch):
                    joint = list(dict.fromkeys(patch + list(EU)))
                    if not all(K in joint for K in missing):
                        full = list(dict.fromkeys(joint + [K for K in K_all if K in patch]))
                        if self.live(full):
                            self.stats['Pr_fail'] += 1; self.failures.append(('Pr', E, EU)); break
                        else:
                            self.stats['joint_dead'] += 1
                    else:
                        self.stats['joint_live'] += 1
            # (Co) coarsening
            for Tp in E:
                if Tp in Sset: continue
                self.stats['Co_checks'] += 1
                for ETp in self.enclosures_given(Tp, patch):
                    joint = list(dict.fromkeys(patch + list(ETp)))
                    csp = self.compat(Tp, set(joint))
                    ok = False
                    if len(csp) == 1:
                        Sp = csp[0][1]
                        ok = self.ch.rel(S, Sp) in self.scaledP
                    if ok:
                        self.stats['joint_live'] += 1; continue
                    # supertiles S (of T) and S' candidates (of T') must be fully present in any tiling
                    # (by (U)+(Pr) applied to T and T'); add their children and re-test consistency
                    if len(csp) == 0:
                        self.stats['joint_dead'] += 1; continue
                    anylive = False
                    for _, Sp in csp:
                        full = list(dict.fromkeys(joint + K_all + self.kids(Sp)))
                        if not self.consistent(full): continue
                        if self.live(full): anylive = True; break
                    if not anylive:
                        self.stats['joint_dead'] += 1; continue
                    if True:
                        self.stats['Co_fail'] += 1; self.failures.append(('Co', E, Tp, ETp, csp)); break
                    else:
                        self.stats['joint_dead'] += 1
            if n % progress == 0:
                log(f"    {n}/{len(encl)} {self.stats} ({time.time()-t0:.0f}s)")
        self.stats['time'] = time.time() - t0
        log("  FINAL:", self.stats)
        return self.stats

def full_run(N, frames, tag, node_limit=20000):
    ch = Chair(N, frames)
    log(f"=== N={N} tag={tag} frames={frames}")
    t0 = time.time()
    P, k = language(ch)
    F = facet_rules(ch, P)
    ntouch, PF = tightness(ch, P, F)
    log(f"  |P|={len(P)} stable at level {k}; |F|={len(F)}; touching poses={ntouch}; |P_F|={len(PF)}; tight={PF == P} ({time.time()-t0:.0f}s)")
    out = dict(N=N, frames=frames, P=sorted(P), F=sorted(F), stable_level=k, touching=ntouch, tight=(PF == P))
    if PF != P:
        out['PF_extra'] = sorted(PF - P); pickle.dump(out, open(f'cert_{tag}.pkl', 'wb'))
        log("  NOT TIGHT -> stop"); return out
    C = Cert(ch, P, node_limit)
    st = C.run()
    out.update(stats=st, failures=C.failures[:50])
    ok = st['U_fail'] == 0 and st['Pr_fail'] == 0 and st['Co_fail'] == 0 and st['limit'] == 0
    out['certified'] = ok
    pickle.dump(out, open(f'cert_{tag}.pkl', 'wb'))
    log(f"=== {tag}: CERTIFIED={ok}  total {time.time()-t0:.0f}s")
    return out

if __name__ == "__main__":
    N = int(sys.argv[1]); frames = eval(sys.argv[2]); tag = sys.argv[3]
    nl = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
    full_run(N, frames, tag, nl)
