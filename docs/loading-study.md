# Phantom loading study: plan, cost and actuals

Simulation only. Coil: the tuned 16-rung low-pass head coil (D 28 cm, L 25 cm, 34 cm closed shield, C = 1.828 pF per
rung, k = 1 at 127.73 MHz unloaded on the production full model). The coil is **not retuned** for the load, so every
case reports its frequency shift. Phantom: a homogeneous cylinder on the coil axis, centred at z = 0, length 20 cm.

## Model

* **Eighth model** (x, y, z ≥ 0), production mesh density (5 / 2.5 / 35 mm, 20 mm inside the phantom), order 2,
  0.159 M unknowns. A centred cylinder has the same three mirror planes as the coil, so the k = 1 member keeps its
  parities: x = 0 PEC, y = 0 PMC, z = 0 PEC (B1 along y), and the partner x = 0 PMC, y = 0 PEC, z = 0 PEC (B1 along
  x). The two runs together give the quadrature B1+ map, as in the symmetry work (`results/lp16_symmetry.json`).
* One mesh per phantom radius. The unloaded reference is the same mesh with the phantom set to εr = 1, σ = 0, so the
  frequency shift is free of mesh-to-mesh noise.
* Three Palace eigen solves per case:
  * `em`: lossless coil (PEC, ideal C), phantom with its conductivity. Its Q is the phantom-only Q.
  * `me`: the partner, same materials; for the quadrature B1+ and E maps.
  * `full`: `em` plus the coil losses: copper (σ = 5.8e7 S/m) on the strips and the shield as a Palace
    `Conductivity` boundary (surface impedance, frequency dependent, so Palace uses its nonlinear eigensolver), and
    a 0.68 Ω series ESR per capacitor, entered as the equivalent parallel R at f0. Its Q is the loaded Q, Q_L. The
    εr = 1, σ = 0 case gives the unloaded Q, Q_U. The strips are internal sheets and Palace applies one surface
    impedance to the sheet, which puts the whole strip current on one face: an upper bound on the strip loss (at most
    2x), so Q_U is on the low side.
* Q = Re f / (2 Im f) from Palace's complex eigenfrequency.
* B1+ normalised to a mean |B1+| of 1 µT over the ROI (r < 8.4 cm, |z| < 4.17 cm), and |E| and the SAR-like
  σ|E|²/2 (W/m³, no tissue density) of the same quadrature field on the z = 0 and y = 0 probe planes (4 mm grid).

## Sweeps

| sweep | points | fixed |
|---|---|---|
| (a) conductivity | σ = 0, 0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0 S/m (12) | εr 78, r 8 cm |
| (b) permittivity | εr = 1, 5, 10, 20, 30, 40, 50, 60, 70, 78 (10; 78 shared with (a)) | σ 0, r 8 cm |
| conductivity only | εr 1, σ 0.6 S/m | r 8 cm |
| (c) radius | r = 5, 6.5, 9.5, 11 cm, each saline (εr 78, σ 0.6) and empty (εr 1) | |

30 cases, 90 solves, run in this order (radius last, so it is the part to cut if the time runs over).

## Checks before the sweep (r = 8 cm, εr 78, σ 0.6 S/m)

1. **Loss tangent instead of conductivity** (speed idea): tan δ = σ / (2π f0 ε0 εr) = 1.0825 at f0 = 127.73 MHz,
   same eighth mesh, `em` run, lossless coil.

   | phantom loss as | k = 1 (MHz) | Q | Palace wall time |
   |---|---|---|---|
   | conductivity 0.6 S/m | 126.867 | 66.80 | 105 s |
   | loss tangent 1.0825 | 126.855 | 65.77 | 150 s |
   | difference | −0.009 % | −1.5 % | |

   They agree within 1.5 % in Q (the loss tangent is exact at f0, and the mode sits 0.7 % below f0), but on this
   model the loss tangent was **not faster**: the eighth-model solve with conductivity took 105 s for 4 modes. The
   slow lossy solves seen before were the full model with 12+ modes. So the sweep uses the conductivity (exact at
   every frequency).
2. **Coil-loss calibration** (εr 1, σ 0, eighth `em`): ESR only Q = 1671 (35 s); copper only Q = 1393 (76 s,
   nonlinear solve). Losses in series: 1 / (1/1671 + 1/1393) = 760.
3. **Full loss in one solve** (saline): Q_L = 60.89 at 126.814 MHz (206 s; one mode at a 0.9 f0 target, because with
   3 modes the nonlinear solver stalled on a strongly damped phantom mode and MUMPS ran out of its memory relaxation).
   Adding the separate parts, 1 / (1/66.80 + 1/760) = 61.4, 0.8 % above the direct value.
4. **Symmetry check**: quarter model (x = 0 PEC / y = 0 PMC and the swapped pair, 0.30 M unknowns) against the eighth
   model on the saline case: f, Q and the quadrature B1+ map. Result below.

## Cost projection (written before the sweep)

Timed on this machine (8 threads, 4 MPI ranks, nothing else running), eighth model, 0.159 M unknowns:

| solve | measured | used for the projection |
|---|---|---|
| `em` or `me`, lossy phantom (σ > 0), 4 modes | 105 s | 105 s |
| `em` or `me`, lossless phantom (σ = 0), 4 modes | (unloaded eighth, 10 modes: 41 s) | 40 s |
| `full`, σ > 0 (nonlinear, 1 mode) | 206 s | 210 s |
| `full`, σ = 0 | (copper only: 76 s) | 100 s |

| group | cases | time per case | total |
|---|---|---|---|
| σ > 0 (11 in (a), conductivity only, 4 saline radii) | 16 | 105 + 105 + 210 s = 7.0 min | 1.87 h |
| σ = 0 (σ = 0 of (a), 9 more in (b), 4 empty radii) | 14 | 40 + 40 + 100 s = 3.0 min | 0.70 h |
| symmetry check (quarter `em` + `me`) | 1 | about 10 min | 0.17 h |
| **total** | | | **about 2.7 h** |

This sits inside the 2 to 7 h budget, so no points are cut. If it runs over, the radius sweep goes first.

## Actuals

Checks before the sweep (free machine): loss tangent 150 s, conductivity 105 s, ESR 35 s, copper 76 s, full loss 258 s
(failed: 3 modes, nonlinear solver stalled on a damped phantom mode) + 206 s (1 mode), eighth partner 102 s, quarter
234 s + 224 s. About 23 min.

Symmetry check (saline, [`results/loading_symmetry_check.json`](../results/loading_symmetry_check.json)):

| | eighth model (x2) | quarter model (x2) |
|---|---|---|
| k = 1 (MHz) | 126.867 / 126.879 | 126.886 / 126.877 |
| Q, phantom only | 66.80 / 66.82 | 66.76 / 66.84 |
| B1+ ROI std/mean | 17.96 % | 17.93 % |
| axial B1+ map, RMS difference | 1.0 % of the mean | |
| unknowns, Palace time | 0.159 M, 102-105 s | 0.303 M, 224-234 s |

Sweep: 90 solves, all converged, 16:02 to 19:08 (3.09 h wall; 3.07 h of Palace time). From 16:05 to 18:55 the CPU was
shared with the builder package test and the unshielded air-box runs (3 more ranks), so most solves ran slower than
the free-machine numbers the projection used (projection: 16 lossy cases x 7.0 min + 14 lossless x 3.0 min + checks = 2.7 h):

| solve | n | projected | measured min / mean / max |
|---|---|---|---|
| `em`, `me`, σ > 0 | 32 | 105 s | 98 / 160 / 374 s |
| `full`, σ > 0 | 16 | 210 s | 114 / 185 / 314 s |
| `em`, `me`, σ = 0 | 28 | 40 s | 28 / 37 / 46 s |
| `full`, σ = 0 | 14 | 100 s | 71 / 135 / 185 s |

Total for the study: about 3.5 h (23 min of checks + 3.1 h sweep), inside the 2 to 7 h budget and 0.8 h over the
2.7 h projection, all of it from the shared CPU. Peak RSS 0.56-0.74 GB per rank.

The full-loss solves need one mode at a 0.9 f0 target; the in-phantom peak |E| and σ|E|²/2 use probes more than 3 mm
inside the phantom (the mesh facets the cylinder, and a probe just inside the surface can land in an air element; at
r = 6.5 cm that gave 206 V/m against 51 V/m 3 mm further in).
