# Census: rep-8 polycubes with ≤ 10 cells and their 2-dissections (session 3c)

Code: `rep_census3.py` (enumerate free 3D polycubes, exact-cover test of 2T by copies of T,
dissection count capped, classes modulo Sym(2T)); logs `rep_census3.log` (≤ 8 cells),
`rep_census3_9_10.log`; follow-ups `rc_Lprism3.log`, `rc_oct8.log`, `rc_cube10.log`, `rc_plate25.log`,
`plate_run.log`.

Free polycubes by size: 1, 1, 2, 7, 23, 112, 607, 3811, 25413, 178083 (sizes 1–10).
Rep-8 (2T tileable by copies of T): 1, 2, 5, 2, 5, 2, 27, 4, 3 (sizes 2–10) — 51 shapes in total.

## Shapes with a unique 2-dissection (classes modulo Sym(2T) = 1)

| cells | shape | hierarchy | verdict |
|---|---|---|---|
| n | 1×1×n bar (all n) | periodic | trivial |
| 8 | 2×2×2 cube | periodic | trivial |
| 9 | 1×3×3 plate | periodic | trivial |
| 7 | **chair** (2×2×2 minus corner) | non-periodic | enforceable (Paper 1) |
| 9 | **L-prism of height 3** (2×2×3 minus a 1×1×3 column) | periodic in z (period 3) | weakly aperiodic at best: the two slabs of every supertile carry the same 2D chair pattern, so the tiling is the 2D chair tiling × a z-period 3 — no periodic P-tiling exists (SAT on 6 tori), consistent with the conjecture, but strong aperiodicity is impossible |

**Result: up to 10 cells the chair is the only rep-8 polycube with a unique 2-dissection and a
non-periodic hierarchy.**

## Shapes with few dissection classes (all fail by Theorem A)

| cells | shape | classes | periodic hierarchy | periodic P-tiling (non-hierarchical) |
|---|---|---|---|---|
| 6 | hex222 (P-prism) | 3 | 1 | 3/3 |
| 6 | hex223 | 64 | 0 | 64/64 |
| 8 | `oct8` (2×2×3 box minus two opposite edge columns, S-shaped) | 10 | 4 | 10/10 |
| 10 | `cube10` (2×2×2 cube + 1×2×1 bump) | 3 | 2 | 2/3 (the third has a periodic hierarchy) |
| 10 | `plate25` (1×2×5 plate) | 11 | 5 | 6/6 non-periodic classes (torus 5×8×8); all 6 fail `mcert` with offsets at level 2 |
| 4 | tripod | 52 | 14 | 52/52 |
| 4 | screw | 10 440 | ≥ 5 (sample) | 40/40 (sample) |

Every non-unique class with a non-periodic hierarchy tested so far (3 + 64 + 6 + 1 + 6 + 38 + 35 = 153
classes) admits a periodic non-hierarchical P-tiling, so by Theorem A no ancestry marking of any
depth enforces it; every unique class either is the chair or has a periodic hierarchy.

## Consequences

1. Within hierarchical (ancestry) markings and 2-dissections, the chair is the unique small 3D
   candidate; the next places to look are (a) inflation factor 3 (rep-27; earlier spot checks of
   bigchair19/cube26 found none, a systematic census ≤ 8 cells is cheap), (b) polycubes with 11–14
   cells (enumeration cost grows ~7× per cell; a C implementation or a direct search restricted to
   shapes inside a 3×3×3 box would be the practical route), (c) 4D (the 4D chair is the only unique
   rep-16 tile ≤ 5 cells; a ≤ 7-cell census is feasible).
2. The periodic-SAT filter (Corollary of Theorem A) decides a class in seconds, so any census can be
   coupled with it directly; `mcert.py` is only needed for the survivors.
