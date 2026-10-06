# Rigidity of hierarchical (ancestry) markings — theorem, proof and census

Session 3b. Code: `rigidity.py` (mixed-substitution test), `persat.py` (periodic SAT search,
inflation, de-substitution), `rigidity_census.py` (per-class census); logs `rc_*.log`, pickles
`rigidity_census_*.pkl`, `persat_*.pkl`.

## 1. Setting and definitions

* T ⊂ Z^N a polycube, D a 2-dissection: 2T = ⋃_{i<2^N} (g_i T + t_i).  The **D-hierarchy** is the
  substitution tiling system obtained by applying D at every level; its hull H(D) is the set of tilings
  all of whose patches occur in supertiles.
* A level-m supertile has shape 2^(m-1)T; level-1 tiles are copies of T.  Poses are (g, τ) with
  g in the hyperoctahedral group and τ a lattice translation (PolyTile/τ convention throughout).
* **Depth-k ancestry marking.**  A tile of a hierarchical tiling is marked by its frame g and the word
  (i_1,…,i_k) where i_j is the index of its level-j ancestor inside its level-(j+1) ancestor.  k = 0 is
  the frame-only marking of Paper 1 (the frames are a choice: *frame variants* of D); k = 1 is the
  "frame + parent index" marking of Session 3.
* **Hierarchical contact language** L_k = {(rel(A,B), m_A, m_B) : A, B face-adjacent marked tiles in
  some tiling of H(D)}.  The **depth-k rule** admits a marked tiling iff every face contact lies in L_k.
  The rule **enforces** the hierarchy if every admissible tiling lies in the (marked) hull.
* The depth-k rule is the one certified in Paper 1 (k=0, chair: 30 poses) and by `mcert.py` (k=1).
  Note that L_k for k ≥ 1 is independent of the frame variant up to relabelling (the index already
  identifies the child), so only dissection *classes* matter for k ≥ 1.

## 2. The inflation lemma

**Lemma 1 (inflation).**  Let 𝒯 be a tiling admitted by the depth-k rule.  Define σ(𝒯) by replacing
every marked tile A = (g, τ, m) by the 2^N marked tiles
  σ_i(A) = (g h_i, 2τ + g τ_i, (i, m_1, …, m_k)),   i < 2^N,
where (h_i, τ_i) are the child maps of D.  Then σ(𝒯) is a tiling admitted by the depth-(k+1) rule.
If 𝒯 is periodic with lattice Λ, σ(𝒯) is periodic with lattice 2Λ.

*Proof.*  σ(𝒯) is a tiling because D dissects 2T exactly and 2·cells(A) + {0,1}^N = cells of the
inflated tile.  Face contacts of σ(𝒯) are of two kinds.
(i) Two children σ_i(A), σ_{i'}(A) of the same tile.  Their relative pose is the D-relative pose of
children i, i', and their marks are (i, m), (i', m).  In the fixed-point tiling every ancestry word m of
length k occurs as the mark of some tile, hence of some level-2 supertile's position in the
depth-(k+1) sense; the children of that supertile realise exactly this contact, so it lies in L_{k+1}.
(ii) Children of two different tiles A, B.  Inflation does not create face adjacencies between tiles
that were not face-adjacent: if 2x+e and 2y+e' (e, e' ∈ {0,1}^N) differ by a unit vector then x − y is
0 or a unit vector.  Hence A and B are face-adjacent in 𝒯 and (r, m_A, m_B) ∈ L_k with r = rel(A,B).
By definition of L_k there is a hierarchical tiling 𝒯_h containing tiles A', B' with rel(A',B') = r and
marks m_A, m_B.  σ(𝒯_h) is again hierarchical, A', B' become level-2 supertiles with children marked
(i, m_A), (i', m_B), and the children of σ(A), σ(B) are the image of the children of σ(A'), σ(B') under
the rigid motion taking A' to A (which takes B' to B because the relative poses agree).  So every
contact between children of A and B occurs in σ(𝒯_h) and lies in L_{k+1}.  Periodicity is clear.  ∎

**Lemma 2 (periodic tilings are outside the hull).**  If the D-hierarchy is not periodic and its hull
is minimal (the substitution on marked prototiles is primitive), then H(D) contains no periodic tiling.
*Proof.*  In a minimal hull every tiling has the same orbit closure; the orbit closure of a periodic
tiling consists of its translates only, so all tilings would be periodic.  ∎
(When minimality is not available we use a direct computational witness instead: a tiling is outside
the hull if some tile has no level-m supertile ancestor for some m — the *de-substitution test*
`desub_Z`.)

## 3. The rigidity theorem for ancestry markings

**Theorem A.**  Let T, D be as above and suppose the frame rule L_0 (for some frame variant) admits a
periodic tiling 𝒯 that is not in H(D).  Then for every k ≥ 0 the depth-k rule admits the periodic
tiling σ^k(𝒯), which is not in H(D) (Lemma 2, or the de-substitution witness, which is inherited by
σ^k(𝒯) since σ maps hull tilings to hull tilings and de-substitution of σ^k(𝒯) recovers 𝒯).
Consequently **no ancestry marking of any depth enforces the D-hierarchy.**  For k ≥ 1 the conclusion
holds for every frame variant of D (relabelling).

*Proof.*  Lemma 1 applied k times.  ∎

**Corollary (necessary condition).**  If some depth-k ancestry rule enforces the D-hierarchy, then
already the frame rule L_0 admits no periodic tiling.  This is a finite SAT problem on tori and is the
cheapest test of a candidate rep-tile.

**Remark (what the theorem does *not* say).**  Goodman-Strauss (1998) shows that every good
substitution tiling can be enforced by *some* finite set of marked tiles; Theorem A only concerns the
natural hierarchical markings (frame + ancestry), which are the ones our certificates handle and the
ones realisable by keys-and-locks on the tile faces without wiring information across the tile.

## 4. Census (all dissection classes, representative frames, tori with side ≤ 8)

| tile | classes | periodic hierarchy | periodic L_0-tiling found | provably non-hierarchical | enforceable? |
|---|---|---|---|---|---|
| chair (7 cells) | 1 | 0 | **0** (tori 4·4·7, 4·7·7, 4·4·14, 7·7·7, 4·7·14) | – | yes (Paper 1, `mcert.py`) |
| hex222 (6) | 3 | 1 (class 0: layered prism, period (0,0,2)) | 3 | 3 (de-substitution fails at level 3) | no, for every depth |
| hex223 (6) | 64 | 0 | 64 | 64 (all fail at level 2) | no, for every depth |
| tripod (4) | 52 | 14 | 52 | 44 (the 8 others have periodic hierarchies) | no, for every depth |
| screw (4) | 40 of 10 440 | 5 | 40 | 39 (the other has a periodic hierarchy) | no (sampled) |
| L-prism-3 (9) | 1 | 1 (z-period 3) | 0 | – | unique dissection but periodic hierarchy |
| oct8 (8) | 10 | 4 | 10 | 10 | no |
| cube10 (10) | 3 | 2 | 2 | 2 | no |
| plate25 (10) | 11 | 5 | 6 of 6 non-periodic | 6 | no |

* "periodic hierarchy": the fixed-point tiling itself is invariant under a translation (verified on
  level-6 supertiles) — such a dissection cannot give aperiodicity by any rule.  This explains the two
  tripod classes (c6, c48) that `mcert.py` left "open": their hierarchies are periodic.
* "periodic L_0-tiling": found by SAT (`persat.py`) on a torus of the listed size; period at most
  (4,4,12).  Every flexible class has one; the chair has none on any tested torus.
* "provably non-hierarchical": the de-substitution test fails in Z^3 (so the tiling is outside the hull
  independent of minimality).  Together with Theorem A this proves that the depth-k rule fails for all
  k ≥ 0 for these classes, replacing the earlier *certificate failures* (which only said the method
  cannot prove enforcement) by actual counterexamples.
* Lemma 1 was also verified computationally: inflating the hex222 counterexamples gives L_1- and
  L_2-admissible tilings with periods (8,8,24) and (16,16,48) (`check_admissible`), which fail
  de-substitution — consistent with the fact that the earlier SAT searches for L_1 on tori of side
  ≤ 12 found nothing.

## 5. What the counterexamples look like, and the failed "mixed substitution" idea

* The natural first guess — replace the dissection by an alternative D' at one level ("mixed
  substitution", `rigidity.py`) — does **not** produce admissible tilings: for every hex222 class, every
  D' (other partitions and single-child frame variants), every depth k ≤ 1 and level M ≤ 4, some
  contact inside the D'-dissected supertile is outside L_k.  The hierarchical contact language is
  strictly finer than "children fit together".
* The actual counterexamples are periodic arrangements of *honest* level-1 or level-2 D-supertiles
  whose mutual contacts are all hierarchical but whose global arrangement is not: for hex222 (classes
  1, 2) the periodic tiling de-substitutes to level 2 (it is a tiling by 2T-supertiles) and fails at
  level 3; for tripod/screw/hex223 most fail already at level 2.  This is the "sliding" mechanism:
  a flexible tile admits several ways to assemble a supertile-sized block, and the rule cannot
  distinguish them because the pieces are identical marked tiles in identical contacts.

## 6. Status of the rigidity conjecture

Conjecture (unchanged): closure/enforcement by ancestry markings ⇔ 2T has a unique dissection into
copies of T (up to the symmetries of 2T).  Evidence now: the "⇒" direction holds in all 159 classes
tested, in the strong form of Theorem A (periodic counterexamples exist whenever the dissection is not
unique, including every tile with a non-periodic hierarchy), and "⇐" holds for the only unique
dissection we know in 3D (the chair) and its N-dimensional analogues (Paper 1).  A proof of "⇒" in
general would need to construct a periodic P-tiling from a second dissection; the mixed-substitution
experiments show this is not a one-level replacement, so a different construction is required.

**Consequence for the programme.**  Within the class of hierarchical (ancestry) markings, strongly
aperiodic rep-2^N polycube monotiles must come from rep-tiles with a *unique* 2-dissection.  The
searches so far (boxes-minus-boxes, all 4D polycubes ≤ 5 cells, 3D tiles ≤ 6 cells tried here) found
no such tile other than the chair and its affine stretches.  The two ways forward are (a) a systematic
census of rep-2^N polycubes with unique dissections (larger cell counts, non-box shapes), using the
periodic-SAT test of the Corollary as the cheap filter, or (b) leaving ancestry markings for
Goodman-Strauss-type markings, which are known to work in principle but are far from a small
keys-and-locks polycube.
