# Solve plan

Machine: 8 threads, 32 GB (WSL sees 25 GB). Another Palace job may share the machine, so at most 4 MPI ranks.

| case | mesh (h_coil / h_gap / h_max, mm) | tets | order | unknowns (est.) | RAM (est.) | wall (est.) | purpose |
|---|---|---|---|---|---|---|---|
| test_coarse | 8 / 4 / 40 | 72 k | 1 | 92 k (actual) | 1.4 GB (actual) | 3 min (actual) | pipeline smoke test, mode ID check |
| lp16_m1p2 | 8 / 4 / 40 | 72 k | 2 | ~0.5 M | 4-6 GB | 10-15 min | p-refinement step |
| lp16_m2 | 5 / 2.5 / 35 | 169 k | 2 | ~1.1 M | 10-15 GB | 25-45 min | h-refinement step; production mesh if converged |
| tune iterations | production mesh | | 2 | same | same | same per iteration | C rescaled by (f_sim/f_target)^2, mesh reused |
| loaded (phantom) | production + 20 mm in phantom | ~250 k | 2 | ~1.6 M | 15-20 GB | 45-90 min | frequency shift, Q |

Unknowns: whether the lumped-C eigenproblem converges with 20 modes from a target of 0.55 f0 (yes on test_coarse: 20 modes);
how far the ladder model is from full wave (test_coarse: k=1 at 112.9 MHz vs 127.7 MHz ladder, -11.6 %);
mesh convergence of the k=1 frequency (decides the production mesh).

Rule: run the cheaper case first; do not start a > 20 min case before the previous level finished and looked sane.
