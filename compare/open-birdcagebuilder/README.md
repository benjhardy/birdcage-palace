# Open BirdcageBuilder vs our ladder model, checked in Palace

Open BirdcageBuilder (D. Brantner, G. Carluccio, C.-L. Chin, C. M. Collins, NYU Langone CAI2R/CBI; ISMRM 2024,
poster #1569; https://github.com/dbnyu/OpenBirdcageBuilder) is the open-source web version of BirdcageBuilder
(Chin, Collins, Li, Dardzinski, Smith, Concepts Magn Reson 2002;15(2):156-163). Both it and our ladder model
(`birdcage_palace/ladder.py`) give a quasi-static first guess for the birdcage capacitors. Here both guesses go into
the same Palace eigenmode model, to see where the k=1 mode lands against the 3 T target.

Simulation only: Palace FEM, unloaded coil, lossless conductors.

## Coil

16 rungs, diameter 28 cm, length 25 cm between end-ring centrelines, rungs 12 mm and end rings 15 mm wide (flat
copper strips, zero thickness), 34 cm closed PEC shield (45 cm long, PEC end caps), 3 T: target 127.7324 MHz
(42.577478 MHz/T x 3 T). Band-pass with C_rung = C_ring (ratio 1.0).

## Open BirdcageBuilder version and licence

- Commit `5a7a3942810a1042b78f2db6b02d5c277a78a3d6` (24 May 2024), version 1.0.1.
- Licence GNU AGPL-3.0. None of its code is copied into this repository. `obcb_harness.mjs` loads its unmodified
  source files from a local clone (outside this repository) into a Node `vm` context and makes the same calls as its
  web page (`get_gui_args()` + `calc_capacitor()` in `BBGuiHandler.js`).
- Harness self-check: the repository's own unit test "App Defaults Tests" expects 205.0515 pF; the harness gives
  205.0514 pF.

## Input mapping

| Open BirdcageBuilder input | value | note |
|---|---|---|
| Number of legs | 16 | |
| Frequency | 127.732434 MHz | |
| Coil radius | 14 cm | rung centre radius |
| Shield radius | 17 cm (shielded) / 0 (no shield) | see "shield" below |
| Leg length | 25 cm | end-ring centreline to centreline, the same length our ladder uses. The clear rung length between the ring inner edges (23.5 cm) gives 1.7 % more C (results JSON, `obcb_leg_length_23p5cm_pf`) |
| Leg type | rectangular, width 1.2 cm | |
| End ring type | rectangular, width 1.5 cm | segment arc length is 2 pi R / N, set by the tool |
| Band-pass | C_leg = C_ring | the tool takes one fixed capacitor and returns the other; the harness finds the value where both are equal by bisection. The "fixed caps on legs" and "fixed caps on ER" modes give the same value (29.977 pF) |

Mismatches that cannot be mapped:

- **Shield.** Open BirdcageBuilder has a shield option: an image rung at radius Rs^2/R, with the same length and the
  opposite current, for every rung. The end rings get no image. Our ladder model images the rungs and the end-ring
  arcs. Neither model has the shield end caps or the finite shield length that Palace has.
- **End-ring coupling.** Open BirdcageBuilder couples the segments of one end ring. Our ladder also couples the two
  end rings to each other.
- **Capacitor gaps.** Neither circuit model has the 5 mm capacitor gaps or the stray capacitance of the strips (to
  each other and to the shield). Palace has both.
- **Conductors.** Both models use a thin-strip self-inductance (Grover; ours via the geometric mean distance 0.2235 w)
  and centre-line filaments for the mutual inductances. Their self-inductances agree: rung 211.5 vs 211.9 nH, ring
  segment 27.4 vs 27.9 nH (unshielded).

## Palace check

Each first-guess C is put into the mesh and Palace configuration of a validated reference case; only the lumped C
values change (`palace_checks.py` asserts that the solver block is identical). Order 2, MUMPS direct solve
("MGMaxLevels": 1), eigen tolerance 1e-8, the same mode count and shift target as the reference runs, 4 MPI ranks,
WSL call in the foreground.

- Low-pass: eighth-symmetry model at production mesh density (`runs/sym_e_m2o2_clean`, 0.156 M unknowns). At the
  Palace-tuned C it gives 127.675 MHz, the same value as the stored reference run (the full model gives 127.734 MHz,
  so the eighth model reads 0.05 % low).
- High-pass and band-pass: the symmetry planes in `birdcage_palace` are implemented for the low-pass coil only, so
  these use the full model on the production mesh (`runs/hp16_3T`, `runs/bp16_3T`, 1.2 M unknowns). The Palace-tuned
  baselines are those reference runs (same mesh, same configuration).

## Results

| | | C (pF) | k=1 in Palace (MHz) | error vs 127.73 MHz | C vs Palace-tuned |
|---|---|---|---|---|---|
| Low-pass | Open BirdcageBuilder | 2.120 | 121.95 | -4.53 % | +16.0 % |
| Low-pass | our ladder first guess | 2.807 | 111.03 | -13.08 % | +53.6 % |
| Low-pass | Palace-tuned | 1.828 | 127.73 | — | — |
| High-pass | Open BirdcageBuilder | 27.857 | 147.79 | +15.70 % | -26.3 % |
| High-pass | our ladder first guess | 36.876 | 129.26 | +1.19 % | -2.4 % |
| High-pass | Palace-tuned | 37.797 | 127.73 | — | — |
| Band-pass | Open BirdcageBuilder | 29.977 | 149.05 | +16.69 % | -27.3 % |
| Band-pass | our ladder first guess | 39.683 | 130.07 (interpolated, see note) | +1.83 % | -3.8 % |
| Band-pass | Palace-tuned | 41.230 | 127.66 | — | — |

Palace k=1 values: low-pass rows from the eighth-symmetry model, except the Palace-tuned row (full model, stored
reference; the eighth model gives 127.675 MHz there). Our low-pass ladder guess on the full model gave 111.08 MHz
(`runs/lp16_final/tune_history_part1.json`), against 111.03 MHz on the eighth model. High-pass and band-pass rows: full
model.

Note on the band-pass ladder row: this value was not solved. To keep the Palace time near 1 h it is interpolated
with 1/f^2 linear in C (the rule `birdcage tune` uses) between the two band-pass production-mesh solves (Open
BirdcageBuilder C and Palace-tuned C). The same rule on the high-pass coil, between the same two kinds of point, gives
129.2554 MHz at the ladder C; the direct solve gives 129.2563 MHz (difference -0.0007 %).

## Without a shield (circuit models only)

| | Open BirdcageBuilder | our ladder | difference |
|---|---|---|---|
| Low-pass C (pF) | 1.783 | 1.804 | -1.13 % |
| High-pass C (pF) | 23.43 | 23.69 | -1.13 % |
| Band-pass C (pF) | 25.21 | 25.50 | -1.13 % |

Without a shield the two circuit models agree to 1.1 % for all three types (the same ratio, because both reduce to one
effective inductance for the k=1 mode). The only unshielded Palace result is low-pass on the coarse order-1 mesh
(`results/shield_compare_coarse.json`): tuned C 0.937 pF, about half of both circuit values. That mesh is not the
production mesh (shielded low-pass: 1.86 pF coarse vs 1.83 pF production), so treat it as approximate.

## Why the guesses differ

1. **The shield is the main difference between the two circuit models.** Without a shield they agree to 1.1 %. With
   the 34 cm shield, Open BirdcageBuilder puts image currents on the rungs only, so its end-ring inductance stays at
   the free-space value. Our ladder also images the end-ring arcs, which lowers the ring inductance. The closed can
   in Palace (with end caps) lowers it too. This is consistent with our high-pass and band-pass guesses landing
   within 2 %, while Open BirdcageBuilder's ring capacitors come out 26-27 % low (k=1 about 16 % high).
2. **Stray capacitance sets the low-pass error.** Neither model has the strips' capacitance to each other and to
   the shield. For the low-pass coil, a fit of 1/f^2 = L_eff (C + C_stray) to three Palace solves (main README) gives
   L_eff 511 nH and C_stray 1.21 pF per rung, about two-thirds of the tuned 1.83 pF. Both circuit guesses are
   therefore too large.
3. **Open BirdcageBuilder is closer for the low-pass coil because two errors partly cancel.** Its guess implies an
   effective inductance of 732 nH (1/(w^2 C) at 2.12 pF), against 553 nH from our ladder and 511 nH from the Palace
   fit. The higher inductance (no end-ring image) offsets most of the missing stray capacitance, so k=1 lands at
   -4.5 % instead of our -13.1 %.
4. **Strip vs filament and mutual inductance are not the cause.** Both use a thin-strip self-inductance and filament
   mutual inductances, and their unshielded self-inductances agree within 2 %. Open BirdcageBuilder does not couple
   the two end rings to each other; ours does. With no shield, all these differences together change C by 1.1 %.

Both tools give a first guess only; Palace retuning (`birdcage tune`) then puts k=1 on target.

## Palace runtime

| case | model | unknowns | wall time |
|---|---|---|---|
| lp_tuned (1.828 pF) | eighth symmetry | 155,824 | 49 s |
| lp_obcb (2.120 pF) | eighth symmetry | 155,824 | 49 s |
| lp_ladder (2.807 pF) | eighth symmetry | 155,824 | 41 s |
| hp_obcb (27.857 pF) | full | 1,195,066 | 973 s |
| bp_obcb (29.977 pF) | full | 1,238,220 | 1607 s |
| hp_ladder (36.876 pF) | full | 1,195,066 | 888 s |

Total 60 min on 4 MPI ranks (WSL, Palace in Apptainer). The full-model solves took longer than the reference
runs (755 s high-pass, 942 s band-pass), probably because of load on the Windows host; WSL showed only Palace running.

## Files

- `obcb_harness.mjs`: runs Open BirdcageBuilder from a local clone. `obcb_cases.json` holds the inputs and
  `obcb_raw.json` the outputs (self-check, and all shielded / unshielded and 25 / 23.5 cm leg-length cases).
- `palace_checks.py`: Palace solves at each C (case folders `runs/obcb_*`, not in git). Results: `palace_checks.json`.
- `make_comparison.py`: writes `results/obcb_comparison.json` and the figure
  (`linkedin-content/birdcage-builder/obcb-comparison.png`).

Reproduce:

    git clone https://github.com/dbnyu/OpenBirdcageBuilder <elsewhere>   # commit 5a7a394
    node compare/open-birdcagebuilder/obcb_harness.mjs <elsewhere>/OpenBirdcageBuilder compare/open-birdcagebuilder/obcb_cases.json > compare/open-birdcagebuilder/obcb_raw.json
    python -u compare/open-birdcagebuilder/palace_checks.py lp_tuned lp_obcb lp_ladder hp_obcb bp_obcb hp_ladder
    python compare/open-birdcagebuilder/make_comparison.py

## Credit

Open BirdcageBuilder: D. Brantner, G. Carluccio, C.-L. Chin, C. M. Collins, "Open BirdcageBuilder", ISMRM 2024,
#1569, https://github.com/dbnyu/OpenBirdcageBuilder (AGPL-3.0). BirdcageBuilder: C.-L. Chin, C. M. Collins, S. Li,
B. J. Dardzinski, M. B. Smith, "BirdcageBuilder: design of specified-geometry birdcage coils with desired current
pattern and resonant frequency", Concepts Magn Reson 2002;15(2):156-163, doi:10.1002/cmr.10030. Open
BirdcageBuilder was run unmodified with the closest inputs we could map; the gaps above come from the model
assumptions listed, not from errors in its code.

