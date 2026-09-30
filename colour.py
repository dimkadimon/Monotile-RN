"""
Corner-colouring model of the marked N-chair (Goodman-Strauss style), general N.

Tile C = [0,2]^N \ [1,2]^N in 'index' coordinates.  Corners: outer corners 2v (v in U) and socket (1,...,1).
At each corner the N facet germs carry a bijection axes -> colours:  m(k) : axes -> colours.
We fix m(0) = m(socket) = identity.  A homochiral frame assignment pi_v (frames h_v = D_s pi_v) determines
m(2v) by the rule "the supertile has the same corner colourings as the tile": child v places its local corner 0
at supertile corner 4v with rotation h_v.
Consistency: every two coincident facet germs (same point, axis, quadrant; opposite sides) have equal colours,
throughout the level-k supertile patch.

Germ of tile with frame g (signed perm), translation (cube-centre units as in chairN: cubes at g(2u+1)+2t):
  corner cube  u (local), corner vertex direction sigma in {+-1}^N (outer corner: the vertex of cube u pointing
  away from the tile; socket: vertex (1,..,1) of the tile).
We generate germs directly from the cube geometry: for a tile, a corner point p with local signs; facet germ
(p, axis i, quadrant q over the other axes, side).  Colour: m(k)(i_local) mapped through the frame.
"""
import itertools, sys, time, pickle
from chairN import *
from cert6 import perm_sign

def perm_of(g):
    """absolute axis j = |g|(i_local): returns dict i_local -> j."""
    perm, signs = g
    # apply(g,x)_j = signs[j] * x[perm[j]]  -> local axis perm[j] goes to absolute axis j
    return {perm[j]: j for j in range(len(perm))}

class ColourModel:
    def __init__(self, N):
        self.N = N
        self.U = [u for u in itertools.product((0, 1), repeat=N) if any(c == 0 for c in u)]
        # local corner data of the tile in cube-centre coordinates (cube centres 2u+1, cube half-width 1):
        # outer corner at vertex 4v (=2*(2v)) i.e. point  (4v_i) ; its cone: tile lies in direction -sigma, sigma_i = +1 if v_i=1 else -1
        self.corners = {}   # key -> (point, sigma, kind)
        for v in self.U:
            p = tuple(4 * x for x in v); sigma = tuple(1 if x else -1 for x in v)
            self.corners[('c', v)] = (p, sigma, 'outer')
        self.corners[('s',)] = (tuple([2] * N), tuple([1] * N), 'socket')  # socket vertex at (2,..,2): missing cube [2,4]^N

    def germs(self, pose, m):
        """list of ((point, axis, quadrant, side), colour) for a tile at pose with corner colouring m (dict key -> perm)."""
        g, t = pose; tt = vscale(2, t); pm = perm_of(g); perm, signs = g
        out = []
        for key, (p, sigma, kind) in self.corners.items():
            rho = m[key]
            P = vadd(apply(g, p), tt)
            # absolute outward signs
            Sig = [0] * self.N
            for i_loc in range(self.N):
                j = pm[i_loc]; Sig[j] = sigma[i_loc] * signs[j]
            for i_loc in range(self.N):
                j = pm[i_loc]
                # germ on hyperplane x_j = P_j; for an outer corner the tile occupies direction -Sig in other axes
                # and lies on side -Sig[j] of the hyperplane. For socket: tile occupies +Sig directions (the missing cube is
                # at +Sig) ... handled by 'kind'.
                if kind == 'outer':
                    quad = tuple(-Sig[l] for l in range(self.N) if l != j); side = -Sig[j]
                else:
                    quad = tuple(Sig[l] for l in range(self.N) if l != j); side = -Sig[j]
                out.append(((P, j, quad, side), rho[i_loc]))
        return out

    def check(self, tiles, m):
        """tiles: list of poses; returns number of colour conflicts among coincident germs (and count of coincidences)."""
        table = {}
        conflicts = 0; coincid = 0
        for T in tiles:
            for (P, j, quad, side), col in self.germs(T, m):
                key = (P, j, quad)
                if key in table:
                    s2, c2 = table[key]
                    if s2 != side:
                        coincid += 1
                        if c2 != col: conflicts += 1
                else:
                    table[key] = (side, col)
        return conflicts, coincid

    def colouring_from_frames(self, ch):
        """derive m(2v) from child frames: child v places its local corner 0 at supertile corner 4v (scaled 2*(2v)... )
        Implemented generically: m(0)=m(socket)=id; for each outer corner K of the supertile find the child germ there
        and read the colours -> m(K)."""
        N = self.N
        m = {key: tuple(range(N)) for key in self.corners}   # start with identity everywhere (placeholder)
        # supertile corners in absolute coords: corner point 2*p for tile corner point p
        # collect germs of children with only m(0)/socket known: use children germs at the supertile corners: the child
        # corner sitting there must be its corner 0 (or socket for T0). Verify and read colours.
        derived = {}
        for key, (p, sigma, kind) in self.corners.items():
            Pbig = vscale(2, p)
            found = None
            for (h, tau) in ch.children:
                pm = perm_of(h); perm, signs = h; tt = vscale(2, tau)
                for k2, (p2, s2, kind2) in self.corners.items():
                    if vadd(apply(h, p2), tt) == Pbig:
                        # absolute signs
                        Sig = [0] * N
                        for i_loc in range(N):
                            j = pm[i_loc]; Sig[j] = s2[i_loc] * signs[j]
                        if tuple(Sig) == sigma and kind2 == kind:
                            found = (k2, h)
            if found is None: raise RuntimeError(f"no child corner at supertile corner {key}")
            k2, h = found
            if k2 not in (('c', (0,) * N), ('s',)):
                return None   # child corner at supertile corner is not corner 0/socket -> not the G-S scheme
            pm = perm_of(h)
            rho = [0] * N
            for i_loc in range(N): rho[pm[i_loc]] = i_loc   # colour i_loc (m(0)=id) sits on absolute axis pm[i_loc]
            derived[key] = tuple(rho)
        return derived

def evaluate_frames(N, frames, levels=(1, 2)):
    ch = Chair(N, frames); cm = ColourModel(N)
    m = cm.colouring_from_frames(ch)
    if m is None: return None
    res = []
    for k in levels:
        res.append(cm.check(ch.patch(k), m))
    return m, res

if __name__ == "__main__":
    # sanity on N=3 solution
    fr30 = {'0': (0, 1, 2), (0, 0, 0): (0, 1, 2), (0, 0, 1): (0, 2, 1), (0, 1, 0): (2, 1, 0), (0, 1, 1): (1, 2, 0),
            (1, 0, 0): (1, 0, 2), (1, 0, 1): (0, 1, 2), (1, 1, 0): (2, 0, 1)}
    print("frames30:", evaluate_frames(3, fr30, levels=(1, 2, 3)))
    import pickle, collections
    res = pickle.load(open('scan3.pkl', 'rb'))
    tab = collections.defaultdict(list)
    for fr, p, cand, odd, new in res:
        r = evaluate_frames(3, fr, levels=(1, 2))
        tab[cand].append(tuple(x[0] for x in r[1]) if r else None)
    for cand in sorted(tab): print(cand, collections.Counter(tab[cand]))
