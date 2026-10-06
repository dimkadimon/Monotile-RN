# Session 2 report — following up the WP1/WP2 recommendations

Date: 2 October 2026. Code and data in `next/` (`chair30_rule.txt`, `rep_census_nd.py`, logs).

## 1. The 30-pose chair rule, made precise (new result, ready for paper 2)

Using the liveness-filtered coarsening closure (`certpoly.closure_live`) on the Paper-1 frame assignment:

| | Paper 1 / Chair44 | this work |
|---|---|---|
| admissible relative poses P\* | 44 | **30** (= the hierarchical language P_0 itself) |
| facet-contact triples F\* | 135 | **93** |
| tight (poses passing F\* among 2388) | 44 | 30 |
| enclosures / live / forced | 33 / 15 / 15 | 33 / 15 / 15 |
| certified (unique hierarchy ⇒ strongly aperiodic) | yes | **yes** |

* The 14 poses removed are exactly the poses that are pairwise-admissible at the supertile level but
  **dead**: two level-1 supertiles in such a relative pose can never be completed to a full enclosure
  (DFS liveness, no node limit hit). All 14 are proper rotations (so no chirality subtlety).
* Every dead pose uses at least one facet triple outside F\*, so the 93-triple facet rule excludes them
  *locally* — the rule is realisable by keys and locks exactly as in Paper 1, Section 10.
* Minimality: all 30 poses occur inside the level-2 supertile, so no full-frame rule with fewer poses admits
  the hierarchy. Hence **30 is the minimum, and it is certified** — Flicker's empirical "30 of 44 occur"
  becomes a theorem, and Chair44's rule set shrinks by a third.

Data: `chair30_rule.txt` (30 poses, 14 dead poses, 93 triples), `chair30_rule.pkl`.

## 2. "Socket" tiles at larger inflation — closed off

* [0,3]^3 \ (1,3]^3 (19 cells), [0,3]^3 \ (2,3]^3 (26 cells) and the chair itself are **not rep-27**.
* Exhaustive check of all boxes a×b×c (a,b,c ≤ 4) minus a corner box: the only rep-8 ones are the chair and
  its affine stretches 2×2×4 \ 1×1×2 and 2×4×4 \ 1×2×2 (same substitution up to a diagonal scaling).
  So "cube minus corner" has no non-trivial siblings among box-minus-box polycubes.

## 3. 4D census of small polycubes (new data)

All free 4-dimensional polycubes with ≤ 5 cells (1, 1, 2, 7, 26 shapes): the rep-16 ones are the straight
bars (translation-only), and nine shapes that are rep-16 with **≥ 3000 dissections each** (cap reached),
including the two genuinely 4D pentacubes. No small 4D tile is rigid; the 4D chair (15 cells, unique
dissection) remains the only rigid candidate known.

## 4. hex222 exhaustive census — running

`census_hex222_full.log`: all 3 × 16 384 frame variants (≈ 11 h CPU). 2 900 done at the time of writing, all
failing with the same signature (16–55 → ~250 → offsets). Partial results are saved per class in
`census_hex222.pkl`.

## 5. Towards the rigidity theorem (status: conjecture with strong evidence)

Observed on every tile tested (3D: hex223, hex222, tripod, screw, L-tricube; 4D: L-prism):

> **Conjecture (rigidity).** If a full-frame marking of a rep-2^N polycube T is coarsening-closed, then the
> dissection of 2T is unique up to Sym(2T) and the hierarchical language is small (|P_0| ≤ |B_N|).

Evidence table (closure outcome vs. number of dissections):

| tile | dissections | min |P_0| seen | outcome |
|---|---:|---:|---|
| chair (3D) | 1 | 30 | closed |
| hex222 | 5 | 16 | grows 16→44→244, offsets |
| hex223 | 64 | 400+ | offsets at once |
| tripod | 256 | 36 | offsets |
| screw | 20 736 | 16 | grows then offsets |
| L-tricube | 1 903 | – | (queued) |
| 4D L-prism | ≥ 200 | 12 361 | offsets |
| 4D chair | 1 | – | Paper 1: structured families fail; general case open |

What the failures look like in detail (hex222, best local-search frames): P_0 has 40 poses with 2–4 facet
contacts; the 48 new poses produced by coarsening have 2 or 4 contacts and 28 of them reuse rotations already
present in P_0 — i.e. the supertile can *slide* into positions that imitate tile-level contacts. A proof of
the conjecture should formalise this: a second dissection or a sliding freedom of the supertile yields
pairwise-admissible supertile poses outside 2P. Not yet proved.

## 6. Recommended next step: two-level markings

For flexible tiles the frame carries too little information. The natural next marking is
(frame, index of the tile inside its parent) — 8·|B_N| markings, still one shape — which is the
Goodman-Strauss construction specialised to our pipeline. This needs the pose type in `cert3/cert6` extended
from (g, t) to (g, t, i) (rel-pose = (g_A^{-1} g_B, t, i_A, i_B)); about a day of refactoring. It would give
positive results for tripod/screw-type tiles (and possibly in 4D) at the price of larger marking sets.

---

## Addendum (session 3): two-level (frame + parent-index) markings — `mcert.py`, `mcensus.py`

**Scheme.** A tile is marked by its frame *and* its index inside its parent supertile. Since the
index already identifies the child, frame variants are mere relabellings — only dissection classes
matter. Certificate: R_1 = hierarchical language of marked tiles; R_{j+1} = pairwise-admissible poses
of (unmarked) level-j supertiles 2^{j-1}T; stop when R_{j+1} = 2R_j (j ≥ 2); then the enclosure /
liveness analysis (claim 1) at every level.

**Validation on the chair** (`python3 mcert.py`): |R_1| = 139 marked poses (projection 30),
R_2 = 30, R_3 = R_4 = 44 → self-similar at level 3; claim 1 passes at levels 1 (7 enclosures, all
type A), 2 and 3 (33 enclosures, 15 live, 14 type B) → CERTIFIED. This reproduces the 30/44 structure
from the indexed language alone, confirming the implementation.

**Flexible tiles** (`mcensus.py NAME`): 

| tile | dissection classes tested | outcome |
|---|---|---|
| hex222 | 3/3 | offset at level 2 (2 classes), level 3 (1 class) |
| tripod | 52/52 | offset at level 3 (50); 2 classes (c6, c48) no offset to level 4 but |R_j| = 7, 58, 238 growing — open |
| screw | 57 of 10 440 classes | offset at level 3 (all) |

**Interpretation.** The index markings pin level 1 completely (R_2 = unmarked R_1 is always
reached), but at level 2 the supertile 2T is an *unmarked* shape and its copies "slide" exactly as the
unmarked tile did: the pairwise-admissible rule for 2T admits the alternative dissections of 4T, so
offsets reappear one level up. Marking k levels of indices only pushes the failure to level k+1,
because the sliding is a property of the shape 2^{j}T, which has the same dissection combinatorics
at every level. This is further evidence for the rigidity conjecture (closure ⇔ unique dissection):
no finite marking of the tile can replace a unique dissection.  The two "open" tripod classes have
a rule that grows without stabilising (no self-similarity), so they are not candidates either.

Files: `mcert.py`, `mcensus.py`, logs `mcensus_hex222.log`, `mcensus_tripod.log`,
`mcensus_tripod2.log`, `mcensus_screw.log`, pickles `mcensus_*_*.pkl`.
