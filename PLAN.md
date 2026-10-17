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
