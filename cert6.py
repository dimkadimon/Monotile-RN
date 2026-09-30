"""
cert6.py -- certificate for markings that are invariant under a subgroup H of the coordinate permutations
(H is a group of symmetries of C_N).  A tile's marking is its frame modulo H.  Rules are pose sets that are
H-bi-invariant:  (g,t) allowed  <=>  (h1 g h2, h1 t) allowed.

Pipeline (as in cert5, with H):
  P0  = hierarchical contact language, symmetrised under H
  P*  = closure under coarsening (pairwise-admissible aligned supertile poses, halved, symmetrised)
  F*  = induced facet rules, tightness check
  Claim 1 (unique compatible supertile + presence) over live enclosures, with tiles canonicalised mod H.
"""
import warnings; warnings.filterwarnings('ignore')
import itertools, pickle, time, sys
from chairN import *
from cert3 import Rules, language, facet_rules, log

def perm_sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]: s = -s
    return s

def subgroup(N, name):
    perms = list(itertools.permutations(range(N)))
    if name == 'trivial': sel = [tuple(range(N))]
    elif name == 'A': sel = [p for p in perms if perm_sign(p) == 1]
    elif name == 'S': sel = perms
    elif name == 'Z':   # cyclic shift group
        sel = [tuple((i + k) % N for i in range(N)) for k in range(N)]
    elif name == 'K4':  # Klein four-group in S_4
        sel = [(0, 1, 2, 3), (1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 1, 0)]
    else: raise ValueError(name)
    return [perm_only(p) for p in sel]

def canonical_frames(N, tau=None):
    """T0 -> id; child v -> id if det(D_s)=+1 else transposition tau."""
    if tau is None: tau = tuple([1, 0] + list(range(2, N)))
    ident = tuple(range(N))
    fr = {'0': ident}
    for v in itertools.product((0, 1), repeat=N):
        if all(v): continue
        neg = sum(v)  # number of -1 signs
        fr[v] = ident if neg % 2 == 0 else tau
    return fr

def symmetrize(P, H):
    out = set()
    for (g, t) in P:
        for h1 in H:
            for h2 in H:
                out.add((compose(compose(h1, g), h2), apply(h1, t)))
    return frozenset(out)

class HCert:
    def __init__(self, ch, H, P, node_limit=50000, use_prop=True):
        self.use_prop = use_prop; self.ch = ch; self.N = ch.N; self.H = H; self.P = frozenset(P)
        self.R = Rules(ch, ch.local_cubes, self.P)
        self.node_limit = node_limit
        self.origin = (ch.id, tuple([0] * ch.N))
        self.kids_local = ch.children
        self.sup_for_child = []
        for (h, tau) in ch.children:
            G = inverse(h); self.sup_for_child.append((G, tuple(-x for x in apply(G, tau))))
        self.scaledP = frozenset((g, vscale(2, t)) for (g, t) in self.P)
        self.stats = dict(limit=0)
        self._canon = {}
        # canonical candidate poses relative to a tile (dedupe by frame mod H)
        self.cand_rel = sorted(set(self.canon_pose(p) for p in self.P))
        # candidates covering each local panel
        self.cover = [[] for _ in self.R.panels]
        for rp in self.cand_rel:
            cb = self.R.cubes(rp)
            for i, (c, d) in enumerate(self.R.panels):
                if vadd(c, d) in cb: self.cover[i].append(rp)

    # ---- canonical frames mod H
    def canon_g(self, g):
        r = self._canon.get(g)
        if r is None:
            r = min(compose(g, h) for h in self.H); self._canon[g] = r
        return r
    def canon_pose(self, pose):
        return (self.canon_g(pose[0]), pose[1])
    def same_tile(self, A, B):
        return A[1] == B[1] and self.canon_g(A[0]) == self.canon_g(B[0])

    def absolute(self, A, rp):
        gA, tA = A; g, t = rp
        return self.canon_pose((compose(gA, g), vadd(tA, apply(gA, t))))

    def kids(self, S):
        G, T = S
        return [self.canon_pose((compose(G, h), vadd(T, apply(G, tau)))) for (h, tau) in self.kids_local]

    def supertiles_containing(self, A):
        """all supertile placements (as (child index, S, kids)) in which tile A (frame mod H) is a child;
        deduplicated by the set of kids."""
        seen = {}
        gA, tA = A
        for h in self.H:
            Ah = (compose(gA, h), tA)
            for c in range(len(self.kids_local)):
                G0, T0 = self.sup_for_child[c]
                S = (compose(Ah[0], G0), vadd(Ah[1], apply(Ah[0], T0)))
                K = self.kids(S)
                key = frozenset(K)
                if key not in seen: seen[key] = (c, S, K)
        return list(seen.values())

    def compat(self, A, patch_tiles):
        """supertiles containing A compatible with the patch: each kid equals a patch tile (mod H) or is disjoint
        from all patch cubes and pairwise-admissible with all patch tiles."""
        occ = set()
        for t in patch_tiles: occ |= self.R.cubes(t)
        pset = set(patch_tiles)
        out = []
        for (c, S, K) in self.supertiles_containing(A):
            ok = True
            for kd in K:
                if kd in pset: continue
                if not self.R.cubes(kd).isdisjoint(occ): ok = False; break
                if not all(self.R.pair_ok(kd, t) for t in patch_tiles): ok = False; break
            if ok: out.append((c, S, K))
        return out

    def candidates_for(self, A, c, d):
        k = self.R.local_panel_index(A, c, d)
        return [self.absolute(A, rp) for rp in self.cover[k]]

    def consistent(self, tiles):
        for i in range(len(tiles)):
            for j in range(i + 1, len(tiles)):
                if not self.R.pair_ok(tiles[i], tiles[j]): return False
        return True

    # ---- enclosures of `centre` given fixed tiles
    def enclosures_given(self, centre, fixed):
        R = self.R
        fixed = [t for t in fixed if t != centre]
        cand = [B for B in (self.absolute(centre, rp) for rp in self.cand_rel)
                if B in fixed or all(R.pair_ok(B, f) for f in fixed)]
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

    # ---- liveness: extension with supertile propagation
    def enclosed(self, A, occ):
        ca = self.R.cubes(A)
        for c in ca:
            for d in self.R.dirs:
                x = vadd(c, d)
                if x not in ca and x not in occ: return False
        return True

    def propagate(self, tiles, occ, decoded):
        if not self.use_prop: return True
        changed = True
        while changed:
            changed = False
            for A in list(tiles):
                if A in decoded or not self.enclosed(A, occ): continue
                opts = self.compat(A, tiles)
                if not opts: return False
                if len(opts) == 1:
                    c, S, K = opts[0]; decoded[A] = S
                    pset = set(tiles)
                    for kd in K:
                        if kd not in pset:
                            tiles.append(kd); pset.add(kd)
                            for x in self.R.cubes(kd): occ[x] = kd
                            changed = True
        return True

    def live(self, patch):
        R = self.R
        tiles = list(patch); occ = {}
        for T in tiles:
            for c in R.cubes(T): occ[c] = T
        decoded = {}
        if not self.propagate(tiles, occ, decoded): return False
        self.nodes = 0; self.limit_hit = False
        res = self._dfs(tiles, occ, decoded, list(patch))
        if self.limit_hit: self.stats['limit'] += 1; return True
        return res

    def _dfs(self, tiles, occ, decoded, must):
        self.nodes += 1
        if self.nodes > self.node_limit: self.limit_hit = True; return False
        R = self.R
        best = None; bestc = None
        for A in must:
            ca = R.cubes(A)
            for c in ca:
                for d in R.dirs:
                    x = vadd(c, d)
                    if x in ca or x in occ: continue
                    cands = [B for B in self.candidates_for(A, c, d) if all(R.pair_ok(B, t) for t in tiles)]
                    if not cands: return False
                    if best is None or len(cands) < len(bestc):
                        best, bestc = (A, c, d), cands
                        if len(cands) == 1: break
                if best is not None and len(bestc) == 1: break
            if best is not None and len(bestc) == 1: break
        if best is None: return True
        for B in bestc:
            t2 = list(tiles); o2 = dict(occ); d2 = dict(decoded)
            t2.append(B)
            for x in R.cubes(B): o2[x] = B
            if not self.propagate(t2, o2, d2): continue
            if self._dfs(t2, o2, d2, must): return True
            if self.limit_hit: return False
        return False


def supertile_poses(ch, R):
    """aligned & offset pairwise-admissible supertile poses (children-level rule R)."""
    N = ch.N
    def kids(G, T):
        return [(compose(G, h), vadd(T, apply(G, tau))) for (h, tau) in ch.children]
    S0 = kids(ch.id, tuple([0] * N))
    # generate candidates: some kid b of S must be P-adjacent to some kid a of S0
    sup_for_child = []
    for (h, tau) in ch.children:
        G0 = inverse(h); sup_for_child.append((G0, tuple(-x for x in apply(G0, tau))))
    seen = set(); cand = set()
    for a in S0:
        for p in R.Plist:
            b = R.absolute(a, p)
            for (G0, T0) in sup_for_child:
                S = (compose(b[0], G0), vadd(b[1], apply(b[0], T0)))
                if S in seen: continue
                seen.add(S)
                K = kids(*S); ok = True
                for x in S0:
                    for y in K:
                        if not R.pair_ok(x, y): ok = False; break
                    if not ok: break
                if ok: cand.add(S)
    return cand

def closure(ch, H, maxit=10):
    P0, k = language(ch, verbose=False)
    P = set(symmetrize(P0, H)); hist = [len(P0), len(P)]
    for it in range(maxit):
        R = Rules(ch, ch.local_cubes, P)
        cand = supertile_poses(ch, R)
        odd = [q for q in cand if any(x % 2 for x in q[1])]
        if odd: return None, hist, f'offset({len(odd)})', cand
        new = set(symmetrize([(g, tuple(x // 2 for x in t)) for g, t in cand], H))
        log(f"  closure it{it}: |P|={len(P)} supertile poses={len(cand)} new={len(new - P)}")
        if new <= P:
            return frozenset(P), hist, ('closed' if new == P else 'closed-strict'), cand
        P |= new; hist.append(len(P))
    return None, hist, 'diverged', None

def tightness_H(ch, P, F):
    allp = ch.all_touching_poses()
    PF = set(p for p in allp if (lambda tr: tr and tr <= F)(ch.panel_triples(p)))
    return len(allp), PF

def claim1(C, progress=100):
    T = C.origin; t0 = time.time()
    encl = C.enclosures_given(T, [])
    st = dict(enclosures=len(encl), live=0, dead=0, U_fail=0, Pr_checks=0, Pr_fail=0, Pr_joint_dead=0,
              Pr_joint_live=0, typeA=0, typeB=0)
    fails = []
    log(f"  enclosures: {len(encl)}")
    for n, E in enumerate(encl):
        patch = [T] + list(E)
        if not C.live(patch): st['dead'] += 1; continue
        st['live'] += 1
        cs = C.compat(T, patch)
        if len(cs) != 1:
            st['U_fail'] += 1; fails.append(('U', E, [c for c, _, _ in cs]))
            if len(fails) <= 3: log(f"    U-fail: {len(cs)} compatible supertiles for enclosure #{n}")
            continue
        c, S, K_all = cs[0]
        missing = [K for K in K_all if K not in patch]
        if not missing: st['typeA'] += 1; continue
        st['typeB'] += 1
        U = K_all[0]
        if U not in patch:
            st['Pr_fail'] += 1; fails.append(('Pr-noU', E)); continue
        st['Pr_checks'] += 1
        for EU in C.enclosures_given(U, patch):
            joint = list(dict.fromkeys(patch + list(EU)))
            if all(K in joint for K in missing): st['Pr_joint_live'] += 1; continue
            if C.live(joint):
                st['Pr_fail'] += 1; fails.append(('Pr', E, EU)); break
            st['Pr_joint_dead'] += 1
        if n % progress == 0: log(f"    {n}/{len(encl)} {st} ({time.time()-t0:.0f}s)")
    st['limit'] = C.stats['limit']; st['time'] = time.time() - t0
    log("  Claim1:", st)
    return st, fails

def run(N, Hname, frames, tag, node_limit=50000, use_prop=True):
    ch = Chair(N, frames); H = subgroup(N, Hname); t0 = time.time()
    log(f"=== N={N} H={Hname} (|H|={len(H)}) tag={tag} frames={frames}")
    Pstar, hist, status, cand = closure(ch, H)
    log(f"  closure: {status} sizes={hist}")
    out = dict(N=N, H=Hname, frames=frames, closure_status=status, closure_hist=hist)
    if Pstar is None:
        pickle.dump(out, open(f'cert6_{tag}.pkl', 'wb')); return out
    scaled = set((g, vscale(2, t)) for g, t in Pstar)
    log(f"  supertile poses == 2P*: {cand == scaled}")
    F = facet_rules(ch, Pstar); ntouch, PF = tightness_H(ch, Pstar, F)
    log(f"  |P*|={len(Pstar)} |F*|={len(F)} touching={ntouch} |P_F*|={len(PF)} tight={PF == Pstar} ({time.time()-t0:.0f}s)")
    out.update(Pstar=sorted(Pstar), F=sorted(F), touching=ntouch, tight=(PF == Pstar), coarsen_exact=(cand == scaled))
    if PF != Pstar:
        out['PF_extra'] = sorted(PF - Pstar); pickle.dump(out, open(f'cert6_{tag}.pkl', 'wb')); return out
    C = HCert(ch, H, Pstar, node_limit=node_limit, use_prop=use_prop)
    st, fails = claim1(C)
    ok = st['U_fail'] == 0 and st['Pr_fail'] == 0 and st['limit'] == 0 and cand == scaled
    out.update(claim1=st, fails=fails[:50], certified=ok)
    pickle.dump(out, open(f'cert6_{tag}.pkl', 'wb'))
    log(f"=== {tag}: CERTIFIED={ok}  ({time.time()-t0:.0f}s)")
    return out

if __name__ == "__main__":
    N = int(sys.argv[1]); Hname = sys.argv[2]; tag = sys.argv[3]
    frames = eval(sys.argv[4]) if len(sys.argv) > 4 else canonical_frames(N)
    run(N, Hname, frames, tag, use_prop=(len(sys.argv) < 6 or sys.argv[5] != 'noprop'))
