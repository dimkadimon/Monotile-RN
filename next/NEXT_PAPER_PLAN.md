# Beyond the chair: candidate tiles for strongly aperiodic monotiles in higher dimensions

Research plan for the second paper in the series (follow-up to arXiv:2610.00916).
Prepared 1 October 2026. Everything below that is marked **[probe]** was computed today with the
scripts in `next/` (`rep8_polycubes.py`, `rep_ndim.py`); everything marked **[lit]** is from the
literature; everything marked **[plan]** is proposed work.

---

## 1. Why look beyond the chair

Paper 1 showed that the Chair44 mechanism is *forced* for the chair: among 2187 homochiral frame
assignments of the 3D chair substitution exactly one (up to conjugation) is coarsening-closed, and
for the 4D chair none of the structured families closes. The chair has a *unique* rep-2^N
dissection, so there is no freedom in the substitution itself — all freedom sits in the child
orientations, and in 4D that freedom seems insufficient.

Other rep-tiles behave differently: many admit **hundreds or thousands of inequivalent
dissections**, each giving a different substitution and hence a different candidate matching
rule. The frame-marking pipeline (hierarchical language → coarsening closure → tightness →
two-shell enclosure analysis) is tile-agnostic; only the `Chair` class in `chairN.py` is specific.
Generalising it to an arbitrary polycube with an arbitrary list of dissections turns the whole
zoo of rep-2^N polycubes into test cases, in every dimension.

Key structural observation: **for a tile with trivial symmetry group a dissection already fixes
every child's frame.** The 12^(2^N-1)-fold frame search of Paper 1 collapses to *one certificate
run per dissection*. Asymmetric rep-tiles are therefore the cheapest candidates in 4D.

---

## 2. The landscape of candidate tile types

| # | Family | Dimension(s) | Substitution | Status | Fit with frame-marking |
|---|--------|-------------|--------------|--------|------------------------|
| A | **Small rep-8 polycubes** (tripod tetracube, screw tetracube, L-tricube prism, L/P-tetracube prisms, two hexacubes) | 3, some generalise to N | rep-2^N, many dissections, rotated children required **[probe]** | not studied as aperiodic monotile candidates | **direct** — same lattice, same group B_N, same pipeline |
| B | **L-tromino prism** C_2 × [0,1]^{N−2} | all N ≥ 3 | rep-2^N with genuinely N-dimensional child rotations **[probe N=4: ≥5 dissections]** | prisms of planar rep-tiles | direct; could give a 4D candidate where the chair failed |
| C | **Chair with other dissections / other inflation** (rep-3^N "L" dissections, mixed inflations) | all N | rep-k^N for k ≥ 3 exists for the chair **[lit]** | matching rules unknown | pipeline needs a general k; contact language larger |
| D | **Pinwheel / quaquaversal tiles** (Conway–Radin right prism, Radin's pinwheel triangle) | 2, 3 | 8 children with irrational relative rotations **[lit]** | matching rules exist for pinwheel (Radin 1994, large tile set) | poses are not lattice-registered → needs a *rotational* version of the contact language (finite local complexity up to rotation holds) |
| E | **Danzer's ABCK tetrahedra / Ammann golden rhombohedra** (icosahedral 3D substitutions) | 3 | 4 resp. 2 prototiles **[lit]** | aperiodic *sets*, not monotiles | frame-marking could ask whether one marked tile suffices via metatile clusters, as the hat does in 2D — speculative |
| F | **Socolar–Taylor and Schmitt–Conway–Danzer** | 3 | SCD: no substitution; ST: half-hex substitution **[lit]** | only *weakly* aperiodic in 3D | not candidates for strong aperiodicity; useful as contrast |
| G | **Translational monotiles à la Greenfeld–Tao** | large d | none (Sudoku-type encoding) **[lit]** | disconnected tiles | different mechanism, out of scope |

Recommendation: the next paper should be about **A + B** (lattice polycube rep-tiles), with C as
an optional extension and D as a stated future direction. That keeps the whole paper within the
certified, computer-checkable framework of Paper 1.

---

## 3. What the probe found today

Enumeration of all free polycubes with ≤ 7 cells (1, 1, 2, 7, 23, 112, 607 shapes) and exact-cover
test of the 2× scaled copy by 8 congruent copies (group B_3, reflections allowed):

| cells | tile | genuinely 3D? | rep-8 dissections | translation-only | needs rotations |
|------:|------|:-:|------:|------:|:-:|
| 1 | cube | – | 1 | 1 | no |
| 2 | domino | prism | 121 | 1 | partly |
| 3 | straight tricube | prism | 1 | 1 | no |
| 3 | **L-tricube** | prism of L-tromino | **1903** | 0 | **yes** |
| 4 | straight tetracube | prism | 1 | 1 | no |
| 4 | L-tetracube | prism | 166 | 0 | yes |
| 4 | square tetracube | prism | 165 | 1 | partly |
| 4 | **tripod tetracube** {0, e1, e2, e3} | **yes** | **256** | 0 | **yes** |
| 4 | **screw tetracube** (chiral) | **yes** | **20736** | 0 | **yes** |
| 5 | straight | prism | 1 | 1 | no |
| 5 | P-pentacube | prism | 82 | 0 | yes |
| 6 | straight | prism | 1 | 1 | no |
| 6 | L-hexacube 2×4 | prism | 135 | 0 | yes |
| 6 | 2×3 block | prism | 153 | 1 | partly |
| 6 | **2×2×3 minus two cells** (no symmetry) | **yes** | **64** | 0 | **yes** |
| 6 | **2×2×2 minus two adjacent cells** | **yes** | **5** | 0 | **yes** |
| 7 | straight | prism | 1 | 1 | no |
| 7 | **chair** (Chair44 shape) | yes | **1** | 0 | yes |

(No other polycube with ≤ 7 cells is rep-8.)

Higher-dimensional probes **[probe]**:

* 4D tripod {0, e1, …, e4}: **not** rep-16 (no dissection).
* 4D L-tromino prism C_2 × [0,1]^2: rep-16, ≥ 5 dissections found immediately, children use
  non-trivial 4D permutations (0,2,1,3), (0,2,3,1), (2,0,3,1) — a genuinely four-dimensional
  substitution with lots of freedom.
* 4D chair: unique dissection (confirms Paper 1).

Reading of the table: the chair is the *most rigid* member of the zoo (1 dissection); the
tripod, screw and the two hexacubes are *flexible* ones (5 … 20736 dissections). The screw tetracube
is chiral, so homochirality — which Paper 1 had to impose — comes for free in any dissection using
one handedness.

---

## 4. Proposed paper 2: "Frame-marking certificates for polycube rep-tiles"

### 4.1 Research questions
1. **Which rep-8 polycubes admit a coarsening-closed frame marking?** I.e. which of the shapes in the
   table are, with suitable facet markings, strongly aperiodic monotiles of R^3? Candidate for the
   headline: *a strongly aperiodic marked tetracube* (4 cells, vs. 7 for Chair44).
2. **Is there a strongly aperiodic marked monotile in R^4?** The L-prism C_2 × [0,1]^2 and other
   asymmetric 4D rep-16 polycubes have many dissections and little frame freedom, so each is a cheap
   certificate run. One success settles the 4D question that Paper 1 left open (for a different tile).
3. **Which dissections force strong aperiodicity and which only weak?** Prisms can give tilings that
   are periodic along the prism axis; the certificate (presence of screw/translation symmetries) tells
   them apart automatically.
4. **Statistics of the hierarchy** (contact counts, enclosures, |Aut| bound) for every success, in the
   format of Table 1 of Paper 1 — a comparative "Chair44 table" for every new tile.

### 4.2 Work packages
| WP | Task | Effort | Tooling |
|----|------|--------|---------|
| 1 | Generalise `Chair` → `PolycubeTile(cells, dissection)`: children given by (g, t) from the exact cover; `cubes`, `contacts`, `panel_triples`, `all_touching_poses` already work on cube sets | 1–2 days | `chairN.py`, `cert6.py` |
| 2 | Dissection enumeration up to symmetry of the scaled tile (quotient the exact covers by Aut(2T) and by Aut(T) on children) | 1 day | `rep8_polycubes.py` |
| 3 | Run closure + tightness + enclosure on every (tile, dissection) for the 3D shapes with ≤ 6 cells; for symmetric tiles add the frame search of Paper 1 (|Aut(T)|^(7) assignments) | hours of CPU | `cert6.py`, `scan3.py` |
| 4 | 4D: enumerate rep-16 polycubes with ≤ 5 cells (small set), all dissections, run certificates | CPU-bound; prune by requiring closure on 3-faces first | new `rep_ndim.py` + WP1 |
| 5 | For every success: facet rule table, keys-and-locks realisation, figures (exploded supertile as in Paper 1 Fig. 2) | 1 day each | `mkfigs.py` |
| 6 | Rigidity results à la Theorem B: among all dissections/frames, how many close? | from WP3 output | — |
| 7 | Optional: rep-27 dissections of the chair and tripod (k = 3 inflation) | later | general k in WP1 |

### 4.3 Expected outcomes / risks
* Best case: a marked tetracube or hexacube strongly aperiodic in R^3 (smaller than Chair44), and a
  strongly aperiodic marked monotile in R^4 from the L-prism family — both new.
* Likely partial case: several 3D successes, 4D still open but with a complete census of small tiles.
* Risk: prisms yield only weak aperiodicity (screw periodicity along the axis). Mitigation: the
  certificate detects this; genuinely 3D shapes (tripod, screw, hexacubes) are not prisms.
* Risk: flexible tiles have huge contact languages (thousands of poses) → closure iteration slower.
  Paper 1's code handled 2388 poses in 2 s; 10^4–10^5 should still be minutes.

### 4.4 Title / framing options
* "Frame-marking certificates for polycube rep-tiles: a census of strongly aperiodic marked
  monotiles in R^3 and R^4"
* "Smaller than the chair: strongly aperiodic marked tetracubes" (if WP3 succeeds for a 4-cell tile)

---

## 5. Immediate next steps (if you say go)
1. WP1 + WP2 (generalised tile class, dissection quotient).
2. First certificate runs on the tripod (256 dissections) and the two hexacubes (64 + 5) — these have
   small or trivial symmetry groups, so the runs are cheap and the answer is decisive.
3. Report back with a table in the style of Paper 1's Table 1 for each success.
