# Solve plan and actuals

Machine: 8 threads, 32 GB (WSL sees 25 GB). Another Palace job sometimes shared the machine, so at most 4 MPI ranks then.

| case | mesh (h_coil / h_gap / h_max, mm) | order | unknowns | estimate | actual | purpose |
|---|---|---|---|---|---|---|
| test_coarse | 8 / 4 / 40 | 1 | 92 k | - | 3.0 min, 1.4 GB, 4 ranks | smoke test, mode ID |
| lp16_m1p2 | 8 / 4 / 40 | 2 | 484 k | 10-15 min, 4-6 GB | 6.6 min, 6.6 GB, 4 ranks (direct MUMPS) | p-refinement |
| lp16_m2 | 5 / 2.5 / 35 | 2 | 1.12 M | 25-45 min, 10-15 GB | 20.4 min, 16.2 GB, 4 ranks | h-refinement, production mesh |
| lp16_m1p2_c2 | 8 / 4 / 40 | 2 | 484 k | 7 min | 6.4 min | 1/f^2 slope vs C for the tune step |
| lp16_final (2 tune solves) | production | 2 | 1.12 M | 20 min each | 14.2 + 13.6 min, 16 GB | k=1 on 127.73 MHz |
| hp16_test | 8 / 4 / 40 | 1 | 100 k | 3 min | 2.2 min, 3 ranks | high-pass functional check |
| lp16_coarse_unloaded / _loaded | 8 / 4 / 40 (+20 mm in phantom) | 2 | 484 k / 498 k | 10-20 min | 3.8 min / 5.4 min (k=1 pair only) | phantom shift and Q |

Lesson: with p-multigrid (default) each shift-and-invert solve at order 2 took about 1 min; `MGMaxLevels = 1`
(direct MUMPS on the full order-2 system) made the whole 20-mode solve 6.6 min. The lossy (phantom) eigenproblem
did not reach 12 modes in 45 min; requesting only the k = 1 pair at tol 1e-6 took 5.4 min.

## Symmetry planes (follow-up)

Same mesh density as the production full model (5 / 2.5 / 35 mm, order 2), C = 1.828 pF, 10 modes per run, target
0.55 f0, MUMPS with `MGMaxLevels = 1`. At most 4 ranks; `pgrep -af palace` checked before every solve.
Coarse smoke tests (8 / 4 / 40 mm, order 1, C = 2.81 pF): quarter 24.8 k unknowns, k=1 at 112.68 MHz
(full 112.85 / 113.04), through-rung quarter 112.78 MHz, eighth 112.52 MHz; 40-75 s each on a shared machine.

| case | model | planes (x / y / z) | tets | unknowns (est.) | estimate | actual | purpose |
|---|---|---|---|---|---|---|---|
| sym_q_em | quarter, planes between rungs | PEC / PMC / - | 44 k | 0.29 M | 2-4 min, 3-4 GB | 71 s, 3.7 GB, 0.30 M (clean) | k=1 member with B1 along y |
| sym_q_me | same mesh | PMC / PEC / - | 44 k | 0.29 M | 2-4 min | 110 s wall (shared) | degenerate partner (B1 along x) |
| sym_q_mm | same mesh | PMC / PMC / - | 44 k | 0.29 M | 2-4 min | 466 s wall (shared with an 8-rank job) | even-k cos family + k=0 |
| sym_q_ee | same mesh | PEC / PEC / - | 44 k | 0.29 M | 2-4 min | 210 s wall (shared) | even-k sin family + k=8 |
| sym_qt_em | quarter, planes through rungs (half rung, C/2) | PEC / PMC / - | 46 k | 0.30 M | 2-4 min | 323 s wall (shared), 0.31 M | validates the C/2 half element |
| sym_qt_em_wrongC | as above, full C on the half element | PEC / PMC / - | 46 k | 0.30 M | 2-4 min | 257 s wall (shared) | shows what a wrong half-element C does |
| sym_e_em | eighth (z=0 PEC, half gap, 2C) | PEC / PMC / PEC | 23 k | 0.15 M | 1-2 min, 2 GB | 139 s Palace (shared), 2.1 GB, 0.16 M | 8x model |
| sym_q_em_n20 | quarter, 20 modes | PEC / PMC / - | 44 k | 0.30 M | 2 min | | same mode count as the full run (fair timing) |
| sym_e_m2o3 | eighth, order 3 | PEC / PMC / PEC | 23 k | 0.45 M | 5-10 min, 6-10 GB | | p-refinement |
| sym_e_m3o2 | eighth, 3.5 / 1.75 / 25 mm, order 2 | PEC / PMC / PEC | 56 k | 0.38 M | 3-6 min, 5 GB | | h-refinement |
| sym_e_m4o2 | eighth, 2.5 / 1.25 / 18 mm, order 2 | PEC / PMC / PEC | 136 k | 0.9 M | 10-20 min, 12-14 GB | | h-refinement (finer than any full model that fits in 32 GB) |
| sym_e_m3o3 | eighth, 3.5 / 1.75 / 25 mm, order 3 | PEC / PMC / PEC | 56 k | 1.1 M | 20-40 min, 16-25 GB | | only if the order-3 memory per unknown allows |

Payoff runs wait for a free machine (no other Palace job) for clean timings and RAM headroom.

## Port tuning with circuitRF (phase 2)

Every capacitor gap becomes a resistive lumped port (R_ref 500 ohm, close to the gap reactance, so S is well
conditioned), Palace driven sweep 100-160 MHz, adaptive fast sweep (PROM, tol 1e-5), each port its own excitation.
Symmetry does not reduce a multi-port excitation by itself (driving one gap breaks every mirror symmetry), but it can be
used explicitly: the four PEC/PMC combinations of the eighth model compute the response to the four symmetry-adapted
excitations, and character orthogonality rebuilds the full 16 x 16 matrix exactly (`birdcage_palace/sparams.py`).
The z = 0 plane is PEC in all four (the z-even family has no gap voltage). So the 16-port costs 4 eighth-model runs of
4 excitations each instead of one full-model run of 16 excitations.

Checks done on a coarse mesh (8 / 4 / 40 mm, order 1): the rebuilt column 1 matches a direct full-model driven solve
to |dS| < 4e-3 on every port (mesh noise; |S| up to 0.8); the S-matrix with every gap closed by 1.828 pF resonates at
127.846 / 128.026 MHz against 127.844 / 128.024 MHz from eigen solves on the same eighth mesh (1e-5). The adaptive
sweep needed 4-5 full solves per excitation; the online (ROM) phase post-processes every output frequency (~0.5 s per
point even on the coarse mesh), so the output grid is 2 MHz (the bare-coil S is so smooth that a cubic spline from a
2 MHz grid is accurate to 1e-8).

| case | model | unknowns | estimate | actual | purpose |
|---|---|---|---|---|---|
| pt_e_{em,me,mm,ee} | eighth, production density, order 2, unloaded, 4 ports each | 0.16 M (complex) | 10-20 min each, 3-4 GB | | 16-port of the bare coil |
| ptl_e_{em,me,mm,ee} | same with the saline phantom | 0.16 M | 10-25 min each, 3-4 GB | | 16-port, loaded |
| verification | eighth eigen with the tuned C (unloaded), loaded eigen 2 modes | 0.16 M | 3-10 min | | the tuned C put k=1 on 127.73 MHz |

circuitRF (1.0.9, headless in WSL): SnP 16-port, every gap closed by its capacitor with ESR 0.68 ohm (capacitor Q about
1000 at 128 MHz), drive gaps 0 and 4 (90 degrees apart) fed through a series match capacitor from 50 ohm ports; goals
|S11|, |S22| <= -20 dB and |S21| <= -20 dB at 127.73 MHz (+-2 kHz unloaded, +-100 kHz loaded).
