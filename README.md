# DDA-FFT — light scattering by arbitrary targets

Light scattering by arbitrarily shaped particles, computed with the **Discrete
Dipole Approximation (DDA)** and **Fast Fourier Transform** acceleration of the
dipole–dipole interaction.

The program solves the linear system of dipole equations iteratively and
reports, as a function of scattering angle,

* extinction and absorption cross sections `Cext`, `Cabs`,
* the 4×4 Mueller matrix of the scattering medium,

averaged over a set of scattering planes and, if requested, over many random
orientations or many random realizations of the particle.

Version 6.3.

---

## Authors

| | |
|---|---|
| Original program | **Zubko E.S.** — Version 6.3, Copyright (C) 2001–2012 |
| Optimization, parallelization and refactoring | **Shmirko K.A.** — <kshmirko@dvo.ru> |

The performance work (what was changed, what was measured, and what was tried
and rejected) is documented in **[REPORT.md](REPORT.md)** (in Russian).
Short version: about **3.4× faster** single-threaded and up to **7.5× faster**
on a multi-core machine for a typical Monte-Carlo run, with output that is
bit-identical to the original code on deterministic targets.

---

## Conventions

* Time dependence of the fields is `exp(−iωt)`.
* Refractive index is given as `m = n′ + i·n″` (the imaginary part is the
  absorption, positive for an absorbing material).
* Phase angles and orientation angles are given in **degrees**; the scattering
  size parameter `x` and the refractive index are dimensionless.
* The target is placed on a cubic lattice of `i_max × j_max × l_max` cells
  (`i_max`, `j_max`, `l_max` should be powers of two); each cell is one dipole.

---

## Building

Requires a C++17 compiler: GCC 7+ or clang 5+. No external libraries.

```sh
sh compile          # g++ -O3 -std=c++17 -pthread *.cpp -o dda_fft.out
```

The `compile` script is a plain two-line wrapper; on a system where `g++` is a
different compiler it can be replaced by any equivalent command, e.g.

```sh
g++ -O3 -std=c++17 -pthread *.cpp -o dda_fft.out
```

---

## Running

### Positional arguments

Every argument may be omitted; the program then asks for it interactively, so
an interactive run is simply `./dda_fft.out`. For scripting, the arguments are
positional. Their meaning is:

| # | Argument | Notes |
|---|---|---|
| 1 | `theta_inf` | bottom limit of the phase-angle range, degrees |
| 2 | `theta_up` | upper limit of the phase-angle range, degrees |
| 3 | `theta_step` | phase-angle increment, degrees |
| 4 | `dst` file | binary state file (running averages; see below) |
| 5 | `dat` file | text result file |
| 6 | number of samples | number of realizations / orientations to average over |
| 7 | orientation mode | `0` — fixed orientation, `1` — Monte-Carlo random orientations |
| 8 | *(mode 0)* `psi` | first rotation, degrees |
| 9 | *(mode 0)* `iota` | second rotation, degrees |
| 10 | *(mode 0)* `kappa` | third rotation, degrees |
| 11 | *(mode 0)* number of scattering planes | |
| 8 | *(mode 1)* number of scattering planes | |
| — | target type | see the table below; this shifts the remaining arguments |
| — | target parameters | depend on the target type (e.g. `m′ m″` for a sphere) |
| — | `i_max` `j_max` `l_max` | lattice dimensions, powers of two |
| — | `x` | size parameter of the target |
| — | polarizability model | `1`…`5`, see below |
| — | `Error_max` | convergence threshold that stops the iteration |
| — | `Max_Step` | maximum number of iterations |
| — | FFT method | `0` — Cooley–Tukey (**recommended**), `1` — Rader–Brenner |
| — | solver | `0`…`3`, see below |
| — | `rTau` | only for solver `3`, the weight factor |

The number of phase angles is `(theta_up − theta_inf)/theta_step + 1`.

### Example: a sphere, 64³ lattice

```sh
# phase angles 0..180 in steps of 10 deg, fixed orientation, 1 sample,
# one scattering plane, sphere, m = 1.5 + 0.01i, 64 x 64 x 64, x = 20,
# Clausius-Mossotti, Error_max = 1e-3, 20 steps, Cooley-Tukey, solver 0
./dda_fft.out 0 180 10 out.dst out.dat 1 0 0 0 0 1 1 1.5 0.01 \
              64 64 64 20 1 0.001 20 0 0
```

For a Monte-Carlo average over 100 random orientations, replace the sample
count and the orientation mode, and give the number of scattering planes:

```sh
./dda_fft.out 0 180 10 out.dst out.dat 100 1 1 1 1.5 0.01 \
              64 64 64 20 1 0.001 20 0 0
```

### Output

`out.dat` is a text file:

```
Cabs=…            Cext=…
  0.000   <16 Mueller matrix elements, row by row>
 10.000   <…>
 …
```

The first line holds the cross sections, then one line per phase angle: the
angle followed by the 16 elements of the 4×4 Mueller matrix in row-major order.
The values are averages over the scattering planes and over all samples.

`out.dst` is a binary file holding the number of completed samples and the
running averages of `Cabs`, `Cext` and the Mueller matrix. **It is rewritten
after every sample**, so an interrupted run can be resumed by starting the
program again with the same file: it continues from the sample count stored in
the file instead of starting over. Delete the file to start from scratch.

### Multi-threading

The samples of a run are independent and are processed by a pool of worker
threads. The default is one worker per CPU core, capped by the number of
samples and by the memory the private workspaces need. To set it explicitly:

```sh
DDA_THREADS=4 ./dda_fft.out …
```

**The best number of workers depends on the lattice size**, because the limit
is memory bandwidth, not cores: on an 8-core Apple M1 the measured optimum was
7 workers for a 64³ lattice but only 2 for 128³. One run with 1, 2, 4 and 8
workers is enough to pick the value for a given machine. The chosen number is
printed at start-up as `Worker threads: N`.

With `N = 1` the program prints the per-iteration discrepancy of every sample;
with `N > 1` it prints one summary line per finished sample instead.

---

## Target types

| # | Target |
|---|---|
| 1 | sphere |
| 2 | ellipsoid |
| 3 | parallelepiped |
| 4 | random particle in a spherical matrix |
| 5 | random particle in an ellipsoidal matrix |
| 6 | random anthrax-like particle |
| 7 | random hedgehog-like particle |
| 8 | random dead anthrax-like particle |
| 9 | random particle in a parallelepiped matrix |
| 10 | random nibble particle in an ellipsoidal matrix |
| 11 | particle read from a file |
| 12 | two spheres in contact along the X axis |
| 13 | hierarchical random particle in an ellipsoidal matrix |
| 14 | sphere with arbitrary radius |

Several materials can be mixed in one target: each material is given the
number of dipoles it contributes and its own refractive index, and the
polarizability is computed per material.

## Polarizability models

| # | Model |
|---|---|
| 1 | Clausius–Mossotti |
| 2 | Clausius–Mossotti plus radiative reaction (Draine) |
| 3 | Clausius–Mossotti plus nonradiative reaction (Goedecke & O'Brien) |
| 4 | Isotropized Lattice Dispersion Relation (Draine) |
| 5 | Lattice Dispersion Relation (Draine) — depends on the incident polarization |

## Iterative solvers

| # | Method |
|---|---|
| 0 | Conjugate gradient (Voevodin) |
| 1 | AA\* minimal-iterations method (Voevodin) |
| 2 | BiCGStab |
| 3 | Simple iteration, needs the weight factor `rTau` |

---

## Notes and known issues

* **Use FFT method `0` (Cooley–Tukey).** Method `1` (Rader–Brenner) performs an
  O(N²) transform per line, so it is already slower on a 32³ lattice and the gap
  grows with size. It is kept only for reference.
* The simple iteration solver (`3`) needs a carefully chosen `rTau`; with a
  weight factor that is too large the process diverges. If no iteration improves
  the discrepancy, the program falls back to the iterate of step 0 and prints
  `Result of the last iteration has been replaced with the most optimal!`.
* `Max_Step` bounds the iteration count; a run converges when the discrepancy
  drops below `Error_max`. On hard targets it may require several hundred steps.
* Memory budget: about **0.34 kB per lattice cell for each worker**, plus
  about 0.38 kB per cell for the interaction matrix, which is shared by all
  workers. A single worker on a 128³ lattice therefore needs roughly 1.4 GB.
  See section 3 of [REPORT.md](REPORT.md) for the exact accounting and for the
  largest lattice that fits in a given amount of RAM.

---

## Files

| File | Contents |
|---|---|
| `main.cpp` | input parsing, sample loop, worker pool, output |
| `dda_func.cpp` | DDA solver, convolutions, scalar products, far field |
| `fft_func.cpp` | three-dimensional FFT |
| `particle.cpp` | particle generation |
| `complex.cpp`, `cvector.cpp` | complex arithmetic, 3-vectors and 3×3 matrices |
| `mmuller.cpp` | Jones and Mueller matrices |
| `randomize.cpp` | random number generation |
| `function.hpp`, `cplx.h`, `cvector.hpp`, `complex.hpp`, `mmuller.hpp`, `randomize.hpp` | declarations |
| `compile` | build script |
| `REPORT.md` | optimization report: profiling, changes, measurements, rejected ideas |
