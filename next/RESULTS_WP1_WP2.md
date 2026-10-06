# WP1 + WP2 report: generalised frame-marking pipeline and first census results

Date: 1 October 2026 (session 1 of paper 2). All code in `next/`; all numbers below were produced
by that code in this session (logs: `census_*.log`, `localsearch_*.log`, pickles `census_*.pkl`).

## 1. What was built (WP1, WP2) — done

* **`polytile.py`** — `PolyTile(cells, dissection, scale=k)`: an arbitrary polycube in Z^N with an arbitrary
  rep-k^N dissection, exposing exactly the interface of `chairN.Chair`, so that `cert3.Rules`,
  `cert3.language`, `cert3.facet_rules`, `cert6.HCert`, `cert6.supertile_poses` run unchanged.
  Also: `dissections()` (exact cover of k·T by congruent copies, all orientations in B_N),
  `dissection_classes()` (quotient by Sym(k·T)), `stabiliser()`, `frame_variants()` (all choices of child
  frames modulo Stab(T), child 0 fixed), `homochiral()`.
* **`certpoly.py`** — `run_tile()`: closure (any inflation k) → facet rules → tightness → two-shell enclosure
  analysis with a *generic* presence test (the pivot child U is no longer the chair's central translate;
  every present child is tried). Also `closure_live()`: a **liveness-filtered closure** (see §3).
* **`census3.py`** — census driver: all dissection classes × all/sampled frame variants of a named shape.
* **`localsearch.py`** — coordinate descent on frames minimising coarsening growth.
* **`rep_ndim.py`, `rep8_polycubes.py`** — rep-tile census in Z^N (from the planning step).

**Validation.** On the chair with the Paper-1 frames the generalised pipeline reproduces Theorem A exactly:
closure 30 → 44 → 44, |F*| = 135, 2388 touching poses, tight, 33 enclosures / 15 live / all forced,
CERTIFIED in 2 s (`python3 certpoly.py`).

## 2. Census results so far (3D, rep-8, full-frame marking, inflation 2)

| tile | cells | Stab | dissection classes (total) | frames tested | closure outcome |
|---|---:|---:|---:|---:|---|
| **hex223** 2×2×3 minus two cells | 6 | 1 | 64 (64) | **all 64 (exhaustive)** | all fail: offset supertile poses already at iteration 0 (|P0| = 400–900) |
| **hex222** 2×2×2 minus a domino | 6 | 4 | 3 (5) | 4000 of 49 152 sampled + local search | all fail; class 0 always passes iteration 0 (|P0| ≈ 16–55) but the rule grows 16 → 44 → 244 and offsets appear at iteration 1–2; classes 1, 2 offset at iteration 0. Best local-search score: 448 new poses after 2 iterations (open, not closed) |
| **tripod** {0,e1,e2,e3} | 4 | 6 | 52 (256) | representative frame of every class + local search (2 restarts × 3 sweeps, first 4 classes) | all fail; ≥ 600 offset poses at iteration 1 even at the local optimum |
| **screw** (chiral tetracube) | 4 | 2 | 10 440 (20 736), 1332 homochiral | 4 frames of each homochiral class (4 436 runs so far, running) | all fail at iteration 1 (growth 16–97 → 130–640, then offsets) |
| **chair** | 7 | 6 | 1 (1) | Paper 1 | certified (control) |
| L-tricube prism (3 cells) | 3 | 4 | 491 (1903), 25 homochiral | queued | – |
| 4D L-prism C_2 × [0,1]^2 | 3 | 16 | ≥ 200 found | 2 frames | |P0| = 12 361, 219 535 offset poses: hopelessly flexible |

Reading: every flexible rep-tile tested so far has a hierarchical contact language that is *not*
self-similar under coarsening — the pairwise-admissible supertile poses are far richer than the tile poses
(growth factors 3–10 per level) and odd-offset supertile placements appear. The chair, with its unique
dissection and small language (30), is so far the only polycube with ≤ 7 cells whose full-frame marking
closes. (Exhaustive only for hex223; hex222/tripod/screw are sampled + locally optimised.)

Other facts established on the way: the chair is **not** rep-27 (3·C_3 has no dissection into 27 chairs),
the screw tetracube and hex222 are rep-27; the 4D tripod is not rep-16.

## 3. A new positive result for the chair: a 30-pose rule

`closure_live` keeps a pairwise-admissible supertile pose only if the two-supertile patch is extendable
(tile-level DFS liveness) and re-checks dropped poses against the final rule — sound, because poses that
occur in a tiling are always extendable. For the chair:

* the 14 poses that Paper 1 added in the coarsening step (30 → 44) are all **dead**: no two supertiles in
  those relative poses can be completed to a full enclosure;
* the fixpoint is **P\* = P_0 = 30 poses, F\* = 93 facet-contact triples** (vs 44 / 135), the rule is tight
  (exactly 30 of the 2388 touching poses pass the 93 triples), and the enclosure analysis certifies it
  (33 enclosures, 15 live, each forcing a unique supertile) — `chair_live: CERTIFIED=True`.

So the marked chair admits a strictly smaller certified facet rule than Chair44's 44-contact rule, and
Flicker's observation that "only 30 of the 44 occur" is upgraded to a certificate: the other 14 are not
merely absent from infinite tilings, they are locally dead.  This is a self-contained new result for paper 2.

## 4. Assessment and next steps

1. The brute-force frame census is CPU-bound (0.4–0.8 s per run; 6^7 frames per tripod class, 4^7 per
   hex222 class). Exhausting hex222 (49k runs, ~10 h) is feasible and would make the negative result for
   6-cell tiles complete; tripod (14.5 M) and the prisms need either symmetry reduction by conjugation or
   a theoretical argument.
2. The failure mode is uniform (growth then offsets), which suggests a **theorem rather than more
   computation**: a necessary condition for coarsening closure in terms of the dissection alone (e.g. the
   supertile must have a "socket" that is filled by a translate, as for the chair). Proving such a
   criterion would explain the census and tell us where to look in 4D.
3. Positive directions that remain open: (a) tiles with a chair-like socket but more cells (e.g.
   3×3×3 minus a 2×2×2 corner — rep-27 — and the "double chair" [0,2]^3 minus two corners), (b) the
   30-pose chair rule and its keys-and-locks realisation (smaller than Chair44), (c) *two-level* markings
   (frame + parent index) for flexible tiles — Goodman-Strauss-type rules, more markings but still one shape.
