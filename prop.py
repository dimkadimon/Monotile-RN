"""Extension search with supertile propagation (decoding)."""
import warnings; warnings.filterwarnings('ignore')
from cert4 import *

class PropCert(Cert):
    def enclosed(self, A, occ):
        ca = self.R.cubes(A)
        for c in ca:
            for d in self.R.dirs:
                x = vadd(c, d)
                if x not in ca and x not in occ: return False
        return True

    def propagate(self, tiles, occ, decoded):
        """decode all fully-enclosed, not-yet-decoded tiles; add forced children. Returns False on contradiction."""
        changed = True
        while changed:
            changed = False
            for A in list(tiles):
                if A in decoded or not self.enclosed(A, occ): continue
                patch = set(tiles)
                opts = []
                for c in range(len(self.kids_local)):
                    S = self.sup_containing(A, c)
                    K = self.kids(S)
                    ok = True
                    for kd in K:
                        if kd in patch: continue
                        if not self.R.cubes(kd).isdisjoint(occ.keys()): ok = False; break
                        if not all(self.R.pair_ok(kd, t) for t in tiles): ok = False; break
                    if ok: opts.append((c, S, K))
                if not opts: return False
                if len(opts) == 1:
                    c, S, K = opts[0]
                    decoded[A] = S
                    for kd in K:
                        if kd not in patch:
                            tiles.append(kd); patch.add(kd)
                            for x in self.R.cubes(kd): occ[x] = kd
                            changed = True
                # else: ambiguous -> leave (sound but weaker)
        return True

    def live(self, patch, node_limit=None):
        R = self.R
        nl = node_limit or self.node_limit
        tiles = list(patch); occ = {}
        for T in tiles:
            for c in R.cubes(T): occ[c] = T
        decoded = {}
        if not self.propagate(tiles, occ, decoded): return False
        must = list(patch)
        self.nodes = 0; self.limit_hit = False
        res = self._dfs(tiles, occ, decoded, must, nl)
        if self.limit_hit:
            self.stats['limit'] += 1; return True
        return res

    def _dfs(self, tiles, occ, decoded, must, nl):
        self.nodes += 1
        if self.nodes > nl: self.limit_hit = True; return False
        R = self.R
        best = None; bestc = None
        for A in must:
            ca = R.cubes(A)
            for c in ca:
                for d in R.dirs:
                    x = vadd(c, d)
                    if x in ca or x in occ: continue
                    cands = [B for B in R.candidates_for(A, c, d) if all(R.pair_ok(B, t) for t in tiles)]
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
            if self._dfs(t2, o2, d2, must, nl): return True
            if self.limit_hit: return False
        return False

if __name__ == "__main__":
    import time
    frames={'0': (1, 2, 0), (0, 0, 0): (1, 2, 0), (0, 0, 1): (0, 2, 1), (0, 1, 0): (0, 2, 1), (0, 1, 1): (0, 1, 2), (1, 0, 0): (0, 2, 1), (1, 0, 1): (1, 2, 0), (1, 1, 0): (2, 0, 1)}
    ch=Chair(3,frames); P,k=language(ch,verbose=False); C=PropCert(ch,P,node_limit=200000)
    T=C.origin
    E=C.enclosures_given(T,[])[0]; patch=[T]+list(E)
    c,S=C.compat(T,set(patch))[0]; K_all=C.kids(S)
    Tp=[t for t in E if t not in set(K_all)][0]
    for ETp in C.enclosures_given(Tp,patch):
        joint=list(dict.fromkeys(patch+list(ETp)))
        for cc,Sp in C.compat(Tp,set(joint)):
            full=list(dict.fromkeys(joint+K_all+C.kids(Sp)))
            if not C.consistent(full): print(cc,"inconsistent"); continue
            t0=time.time(); r=C.live(full); print(cc, "in2P",ch.rel(S,Sp) in C.scaledP, "live",r,"limit",C.limit_hit,"nodes",C.nodes,f"{time.time()-t0:.1f}s")
