# DDA-FFT — light scattering by arbitrary targets

**English** | [Русский](#русский)

Light scattering by arbitrarily shaped particles, computed with the **Discrete
Dipole Approximation (DDA)** and **Fast Fourier Transform** acceleration of the
dipole–dipole interaction.

The program solves the linear system of dipole equations iteratively and
reports, as a function of scattering angle,

* extinction and absorption cross sections `Cext`, `Cabs`,
* the 4×4 Mueller matrix of the scattering medium,

averaged over a set of scattering planes and, if requested, over many random
orientations or many random realizations of the particle.

Version 7.2.

---

## Authors

| | |
|---|---|
| Original program | **Zubko E.S.** — Version 6.3, Copyright (C) 2001–2012 |
| Version 7.2 — modification, optimisation and the command-line interface | **Shmirko K.A. (K.A.Shmirko)** — <kshmirko@dvo.ru>, 2026 |

The performance work (what was changed, what was measured, and what was tried
and rejected) is documented in **[REPORT.md](REPORT.md)** (in Russian).
Short version: about **4.75× faster** single-threaded on a 64³ lattice (3.35×
from the first four stages plus 1.42× from an own mixed-radix FFT kernel), and
up to **7.5× faster** on a multi-core machine for a typical Monte-Carlo run
(measured before the kernel rewrite; the kernel adds about the same 1.4× at
four workers). Through the first four stages the output was bit-identical to
the original code; with the new kernel it agrees to `float` precision — about
1e-4 relative on the weakest significant Mueller-matrix elements.

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

Requires a C++17 compiler: GCC 7+ or clang 5+. The command-line parser is the
vendored argtable3 (**BSD-3-Clause**, `third_party/argtable3/`), so there are
no external dependencies to install.

```sh
sh compile
```

The `compile` script builds argtable3 as C and links it into the program:

```sh
cc -O2 -std=c11 -c third_party/argtable3/argtable3.c -o third_party/argtable3/argtable3.o
g++ -O3 -std=c++17 -pthread *.cpp third_party/argtable3/argtable3.o -o dda_fft.out
```

On a system where `g++` is a different compiler the second command can be
replaced by any equivalent one.

### Git hooks

The repository ships a `pre-commit` hook in `.githooks/` that rejects commits
whose author or committer is the agent identity `Agent <agent@local>` instead
of the configured author. Enable it once per clone:

```sh
git config core.hooksPath .githooks
```

To bypass the check for a single commit (not recommended):
`AGENT_AUTHOR_OK=1 git commit ...`.

---

## Running

### Options

The program is configured with named options (parsed with argtable3).
`./dda_fft.out --help` prints the same list. Everything that is missing is
asked for interactively, so a bare `./dda_fft.out` still walks through the
whole questionnaire; `--interactive` forces that when some options are given.

| Option | Meaning | Default |
|---|---|---|
| `-h, --help` | print the options and exit | |
| `-V, --version` | print the version and exit | |
| `-i, --interactive` | ask for every value that is not given | off |
| `-o, --output <file>` | text result file | `out.dat` |
| `-D, --state <file>` | binary state file, resumed if it exists | `out.dst` |
| `--restart` | delete the state file and start from scratch | off |
| `-n, --samples <N>` | number of samples / orientations | `1` |
| `-a, --theta-min <deg>` | bottom limit of the phase angles | `0` |
| `-b, --theta-max <deg>` | upper limit of the phase angles | `180` |
| `-s, --theta-step <deg>` | phase-angle increment | `10` |
| `-O, --orientation <mode>` | `fixed` or `mc` | `fixed` |
| `--psi --iota --kappa <deg>` | rotations, fixed orientation only | `0` |
| `-p, --planes <N>` | number of scattering planes | `1` |
| `-t, --target <name\|1..14>` | target type, see the table below | required |
| `-g, --grid <i> [j l]` | lattice dimensions, powers of two; one value means a cube | required |
| `-x, --size-parameter <x>` | size parameter of the target | required |
| `-m, --material <mre> <mi>` | refractive index of one material, repeat for a mixture | required |
| `--seeds <s>...` | random seeds, one per material | random |
| `--aspect <b/a> <c/a>` | aspect ratios (targets 2, 3, 5, 9, 10, 13) | |
| `--radius-ratio <R1/R2>` | radius ratio (target 12) | |
| `--levels <N>` | number of hierarchical levels (target 13) | |
| `--bites <N>`, `--depth <r>` | surface layer (target 10) | |
| `--packing <rho>` | volume density (targets 6, 7, 8, 14) | |
| `--monomer-radius <r>`, `--monomer-length-ratio <l/r>` | monomer shape (targets 6, 7, 8) | |
| `-f, --shape-file <file>` | file with the target shape (target 11) | |
| `-P, --polarizability <model>` | `1`…`5` or `cm,cmrr,cmnr,ildr,ldr` | `1` |
| `-e, --tolerance <err>` | convergence threshold `Error_max` | `1e-3` |
| `-S, --max-steps <N>` | maximum number of iterations | `20` |
| `-F, --fft <method>` | `0`/`ct` Cooley–Tukey, `1`/`rb` Rader–Brenner | `0` |
| `-M, --solver <method>` | `0`…`3` or `cg,aa,bicgstab,simple` | `0` |
| `--tau <factor>` | weight factor, required by solver `3` | |
| `-j, --threads <N>` | worker threads, `0` = auto | `0` |
| `--seed <u64>` | random seed for reproducible runs | random |

The number of phase angles is `(theta_max − theta_min)/theta_step + 1`.
Grid dimensions must be powers of two; an invalid size is rejected instead of
silently reading past the FFT buffers.

List options accept both forms, so these are equivalent:

```sh
./dda_fft.out --material 1.5 0.01 --grid 64 64 64
./dda_fft.out --material 1.5,0.01 --grid 64,64,64
```

### Example: a sphere, 64³ lattice

```sh
# phase angles 0..180 in steps of 10 deg, fixed orientation, 1 sample,
# one scattering plane, sphere, m = 1.5 + 0.01i, 64 x 64 x 64, x = 20,
# Clausius-Mossotti, Error_max = 1e-3, 20 steps, Cooley-Tukey, solver 0
./dda_fft.out --target sphere --material 1.5 0.01 --grid 64 --size-parameter 20
```

For a Monte-Carlo average over 100 random orientations, ask for the averaging
mode and give the number of scattering planes:

```sh
./dda_fft.out --samples 100 --orientation mc --planes 1 \
              --target sphere --material 1.5 0.01 --grid 64 --size-parameter 20
```

A mixture is given by repeating `--material`; the random targets also need one
`--seeds` value per material:

```sh
./dda_fft.out --target random-sphere --material 1.5 0.01 --material 1.33 0.0 \
              --seeds 1234 5678 --grid 64 --size-parameter 20
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
./dda_fft.out --threads 4 …
```

`DDA_THREADS` in the environment is still honoured when `--threads` is not
given.

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
* The 1-D FFT kernel comes from `fft_kernel.hpp` (mixed-radix 4/8, in place)
  when the line length is a power of two; otherwise `fft_func.cpp` falls back to
  the original radix-2 kernel. Building with `-DFFT_KERNEL_LEGACY` forces the
  original kernel everywhere, which reproduces the pre-rewrite output bit for
  bit.
* The simple iteration solver (`3`) needs a carefully chosen `rTau`; with a
  weight factor that is too large the process diverges. If no iteration improves
  the discrepancy, the program falls back to the iterate of step 0 and prints
  `Result of the last iteration has been replaced with the most optimal!`.
* `Max_Step` bounds the iteration count; a run converges when the discrepancy
  drops below `Error_max`. On hard targets it may require several hundred steps.
* **Particle generation can hang** for the random worm-like targets (6, 7, 8) and
  for the hierarchical target (13) when the requested geometry does not fit: the
  worm placement and the seed-cell search of `particle.cpp` retry without an
  attempt limit. Pick `--packing`, `--monomer-*` and the lattice size so that the
  particle fits, and use an external timeout in scripts. See section 7.6 of
  [REPORT.md](REPORT.md).
* Memory budget: about **0.34 kB per lattice cell for each worker**, plus
  about 0.38 kB per cell for the interaction matrix, which is shared by all
  workers. A single worker on a 128³ lattice therefore needs roughly 1.4 GB.
  See section 3 of [REPORT.md](REPORT.md) for the exact accounting and for the
  largest lattice that fits in a given amount of RAM.

---

## Files

| File | Contents |
|---|---|
| `main.cpp` | sample loop, worker pool, output |
| `options.cpp` | command-line parsing (argtable3), validation, interactive prompts |
| `dda_func.cpp` | DDA solver, convolutions, scalar products, far field |
| `fft_func.cpp` | three-dimensional FFT |
| `fft_kernel.hpp` | mixed-radix 4/8 in-place 1-D FFT kernel used by `fft_func.cpp` |
| `particle.cpp` | particle generation |
| `complex.cpp`, `cvector.cpp` | complex arithmetic, 3-vectors and 3×3 matrices |
| `mmuller.cpp` | Jones and Mueller matrices |
| `randomize.cpp` | random number generation |
| `function.hpp`, `cplx.h`, `cvector.hpp`, `complex.hpp`, `mmuller.hpp`, `randomize.hpp`, `options.hpp` | declarations |
| `compile` | build script |
| `third_party/argtable3/` | vendored argtable3 v3.3.1 (BSD-3-Clause), the CLI parser |
| `bench/` | standalone harness for the 1-D kernel: delta tests, accuracy, A/B timing |
| `REPORT.md` | optimization report: profiling, changes, measurements, rejected ideas |

---

<a id="русский"></a>

# DDA-FFT — рассеяние света на частицах произвольной формы

Рассеяние света частицами произвольной формы, рассчитываемое методом
**дискретных диполей (Discrete Dipole Approximation, DDA)** с ускорением
диполь-дипольного взаимодействия через **быстрое преобразование Фурье (FFT)**.

Программа итерационно решает систему уравнений для дипольных моментов и выдаёт в
зависимости от угла рассеяния:

* сечения экстинкции и поглощения `Cext`, `Cabs`,
* матрицу Мюллера 4×4 рассеивающей среды,

усреднённые по набору плоскостей рассеяния и, при необходимости, по множеству
случайных ориентаций или множеству случайных реализаций частицы.

Версия 7.2.

---

## Авторы

| | |
|---|---|
| Исходная программа | **Зубко Е.С. (Zubko E.S.)** — Version 6.3, Copyright (C) 2001–2012 |
| Версия 7.2 — изменение, оптимизация и командный интерфейс | **Шмирко К.А. (K.A.Shmirko)** — <kshmirko@dvo.ru>, 2026 |

Работа по производительности (что изменено, что измерено и что проверено и
отвергнуто) описана в **[REPORT.md](REPORT.md)**. Кратко: примерно **4.75×**
ускорения в одном потоке на решётке 64³ (3.35× дали первые четыре этапа и ещё
1.42× — собственное смешанное radix-4/8 ядро FFT) и до **7.5×** на
многоядерной машине на типичном прогоне с усреднением по реализациям (замер до
переписывания ядра; ядро добавляет примерно те же 1.4× на четырёх потоках). До
этапа с ядром результат был побитово равен исходному коду; с новым ядром он
совпадает с точностью `float` — около 1e-4 по относительной величине на самых
слабых значимых элементах матрицы Мюллера.

---

## Соглашения

* Временная зависимость полей — `exp(−iωt)`.
* Показатель преломления задаётся как `m = n′ + i·n″` (мнимая часть отвечает за
  поглощение, положительна для поглощающего материала).
* Фазовые углы и углы ориентации задаются в **градусах**; параметр размера `x` и
  показатель преломления безразмерны.
* Частица размещается на кубической решётке `i_max × j_max × l_max` ячеек
  (`i_max`, `j_max`, `l_max` должны быть степенями двойки); каждая ячейка — один
  диполь.

---

## Сборка

Требуется компилятор с C++17: GCC 7+ или clang 5+. Разбор командной строки —
вендоренная argtable3 (**BSD-3-Clause**, `third_party/argtable3/`), так что
внешних зависимостей устанавливать не нужно.

```sh
sh compile
```

Скрипт `compile` собирает argtable3 как C и линкует её с программой:

```sh
cc -O2 -std=c11 -c third_party/argtable3/argtable3.c -o third_party/argtable3/argtable3.o
g++ -O3 -std=c++17 -pthread *.cpp third_party/argtable3/argtable3.o -o dda_fft.out
```

На системе, где `g++` означает другой компилятор, вторую команду можно заменить
эквивалентной.

### Git-хуки

В репозитории есть хук `pre-commit` в каталоге `.githooks/`, который
запрещает коммиты, у которых автор или коммиттер — служебная идентичность
`Agent <agent@local>` вместо настроенного автора. Включите его один раз
после клонирования:

```sh
git config core.hooksPath .githooks
```

Обойти проверку для одного коммита (не рекомендуется):
`AGENT_AUTHOR_OK=1 git commit ...`.

---

## Запуск

### Опции

Программа настраивается именованными опциями (разбор — argtable3). Тот же
список печатает `./dda_fft.out --help`. Всё, что не задано, спрашивается в
интерактивном режиме, поэтому запуск без аргументов по-прежнему проходит весь
опросник; `--interactive` включает его принудительно, если часть опций задана.

| Опция | Смысл | По умолчанию |
|---|---|---|
| `-h, --help` | напечатать опции и выйти | |
| `-V, --version` | напечатать версию и выйти | |
| `-i, --interactive` | спрашивать всё, что не задано | выкл. |
| `-o, --output <файл>` | текстовый файл результата | `out.dat` |
| `-D, --state <файл>` | двоичный файл состояния, продолжается, если существует | `out.dst` |
| `--restart` | удалить файл состояния и начать с нуля | выкл. |
| `-n, --samples <N>` | число образцов / ориентаций | `1` |
| `-a, --theta-min <град>` | нижняя граница фазовых углов | `0` |
| `-b, --theta-max <град>` | верхняя граница фазовых углов | `180` |
| `-s, --theta-step <град>` | шаг по фазовому углу | `10` |
| `-O, --orientation <режим>` | `fixed` или `mc` | `fixed` |
| `--psi --iota --kappa <град>` | повороты, только для фиксированной ориентации | `0` |
| `-p, --planes <N>` | число плоскостей рассеяния | `1` |
| `-t, --target <имя\|1..14>` | тип частицы, см. таблицу ниже | обязательно |
| `-g, --grid <i> [j l]` | размеры решётки, степени двойки; одно значение — куб | обязательно |
| `-x, --size-parameter <x>` | параметр размера частицы | обязательно |
| `-m, --material <mre> <mi>` | показатель преломления материала; повторяется для смеси | обязательно |
| `--seeds <s>...` | случайные seed'ы, по одному на материал | случайно |
| `--aspect <b/a> <c/a>` | отношения осей (типы 2, 3, 5, 9, 10, 13) | |
| `--radius-ratio <R1/R2>` | отношение радиусов (тип 12) | |
| `--levels <N>` | число иерархических уровней (тип 13) | |
| `--bites <N>`, `--depth <r>` | поверхностный слой (тип 10) | |
| `--packing <rho>` | объёмная плотность (типы 6, 7, 8, 14) | |
| `--monomer-radius <r>`, `--monomer-length-ratio <l/r>` | форма мономера (типы 6, 7, 8) | |
| `-f, --shape-file <файл>` | файл с формой частицы (тип 11) | |
| `-P, --polarizability <модель>` | `1`…`5` или `cm,cmrr,cmnr,ildr,ldr` | `1` |
| `-e, --tolerance <err>` | порог невязки `Error_max` | `1e-3` |
| `-S, --max-steps <N>` | максимальное число итераций | `20` |
| `-F, --fft <метод>` | `0`/`ct` Кули–Тьюки, `1`/`rb` Радер–Бреннер | `0` |
| `-M, --solver <метод>` | `0`…`3` или `cg,aa,bicgstab,simple` | `0` |
| `--tau <множитель>` | весовой множитель, нужен методу `3` | |
| `-j, --threads <N>` | рабочих потоков, `0` = авто | `0` |
| `--seed <u64>` | seed для воспроизводимых прогонов | случайно |

Число фазовых углов равно `(theta_max − theta_min)/theta_step + 1`.
Размеры решётки обязаны быть степенями двойки; неверный размер теперь
отвергается, а не приводит к чтению за границей FFT-буферов.

Списочные опции принимают обе формы записи, они эквивалентны:

```sh
./dda_fft.out --material 1.5 0.01 --grid 64 64 64
./dda_fft.out --material 1.5,0.01 --grid 64,64,64
```

### Пример: сфера, решётка 64³

```sh
# фазовые углы 0..180 с шагом 10°, фиксированная ориентация, 1 образец,
# одна плоскость рассеяния, сфера, m = 1.5 + 0.01i, 64 x 64 x 64, x = 20,
# Клаузиус–Моссотти, Error_max = 1e-3, 20 шагов, Кули–Тьюки, метод решения 0
./dda_fft.out --target sphere --material 1.5 0.01 --grid 64 --size-parameter 20
```

Для усреднения по 100 случайным ориентациям нужно задать режим усреднения и
число плоскостей рассеяния:

```sh
./dda_fft.out --samples 100 --orientation mc --planes 1 \
              --target sphere --material 1.5 0.01 --grid 64 --size-parameter 20
```

Смесь задаётся повторением `--material`; случайным типам нужен ещё один seed
на каждый материал:

```sh
./dda_fft.out --target random-sphere --material 1.5 0.01 --material 1.33 0.0 \
              --seeds 1234 5678 --grid 64 --size-parameter 20
```

### Вывод

`out.dat` — текстовый файл:

```
Cabs=…            Cext=…
  0.000   <16 элементов матрицы Мюллера, построчно>
 10.000   <…>
 …
```

Первая строка — сечения, далее по строке на каждый фазовый угол: угол и 16
элементов матрицы Мюллера 4×4 в порядке строка за строкой. Значения усреднены по
плоскостям рассеяния и по всем образцам.

`out.dst` — двоичный файл, в котором хранятся число завершённых образцов и текущие
средние `Cabs`, `Cext` и матрицы Мюллера. **Он перезаписывается после каждого
образца**, поэтому прерванный прогон можно продолжить, запустив программу снова с
тем же файлом: она продолжит со счётчика образцов, записанного в файл, а не
начнёт заново.
Чтобы начать с нуля, файл нужно удалить.

### Многопоточность

Образцы в прогоне независимы и обрабатываются пулом рабочих потоков. По умолчанию
берётся по одному потоку на ядро, но не больше числа образцов и не больше, чем
позволяет память под приватные рабочие массивы. Явно:

```sh
./dda_fft.out --threads 4 …
```

Переменная окружения `DDA_THREADS` по-прежнему учитывается, если `--threads`
не задан.

**Лучшее число потоков зависит от размера решётки**, потому что ограничение —
полоса пропускания памяти, а не ядра: на 8-ядерном Apple M1 измеренный оптимум
составил 7 потоков для 64³, но всего 2 для 128³. Одного прогона с 1, 2, 4 и 8
потоками достаточно, чтобы подобрать значение под конкретную машину. Выбранное
число печатается при старте строкой `Worker threads: N`.

При `N = 1` программа печатает невязку на каждой итерации каждого образца; при
`N > 1` — одну сводку на завершённый образец.

---

## Типы частиц

| № | Частица |
|---|---|
| 1 | сфера |
| 2 | эллипсоид |
| 3 | параллелепипед |
| 4 | случайная частица в сферической матрице |
| 5 | случайная частица в эллипсоидальной матрице |
| 6 | случайная частица типа «антракс» |
| 7 | случайная частица типа «ёж» |
| 8 | случайная «мёртвая» частица типа «антракс» |
| 9 | случайная частица в параллелепипедной матрице |
| 10 | случайная «изъеденная» частица в эллипсоидальной матрице |
| 11 | частица, прочитанная из файла |
| 12 | две соприкасающиеся сферы (вдоль оси X) |
| 13 | иерархическая случайная частица в эллипсоидальной матрице |
| 14 | сфера произвольного радиуса |

В одной частице можно смешивать несколько материалов: для каждого задаётся число
составляющих его диполей и свой показатель преломления, поляризуемость считается
для каждого материала отдельно.

## Модели поляризуемости

| № | Модель |
|---|---|
| 1 | Клаузиуса–Моссотти |
| 2 | Клаузиуса–Моссотти плюс радиационная реакция (Draine) |
| 3 | Клаузиуса–Моссотти плюс нерадиационная реакция (Goedecke & O'Brien) |
| 4 | изотропизированное соотношение дисперсии решётки (Draine) |
| 5 | соотношение дисперсии решётки (Draine) — зависит от поляризации падающего поля |

## Итерационные методы

| № | Метод |
|---|---|
| 0 | сопряжённых градиентов (Воеводин) |
| 1 | метод минимальных итераций AA\* (Воеводин) |
| 2 | BiCGStab |
| 3 | простая итерация, требует весового множителя `rTau` |

---

## Замечания и известные особенности

* **Используйте метод FFT `0` (Кули–Тьюки).** Метод `1` (Радера–Бреннера)
  выполняет O(N²) преобразование на строку, поэтому он медленнее уже на решётке
  32³, и разрыв растёт с размером. Он оставлен только для справки.
* 1-D ядро БПФ берётся из `fft_kernel.hpp` (смешанное radix-4/8, на месте), когда
  длина строки — степень двойки; иначе `fft_func.cpp` возвращается к исходному
  radix-2 ядру. Сборка с `-DFFT_KERNEL_LEGACY` принудительно использует исходное
  ядро — оно воспроизводит результат до переписывания побитово.
* Метод простой итерации (`3`) требует аккуратно подобранного `rTau`; при слишком
  большом весовом множителе процесс расходится. Если ни одна итерация не уменьшила
  невязку, программа возвращается к итерации шага 0 и печатает
  `Result of the last iteration has been replaced with the most optimal!`.
* `Max_Step` ограничивает число итераций; прогон считается сошедшимся, когда
  невязка опускается ниже `Error_max`. На трудных частицах может потребоваться
  несколько сотен шагов.
* **Генерация частицы может зависнуть** для случайных червеобразных целей (6, 7, 8) и
  для иерархической (13), когда запрошенная геометрия не помещается: размещение
  червей и поиск затравок в `particle.cpp` повторяются без ограничения числа
  попыток. Подбирайте `--packing`, `--monomer-*` и размер решётки так, чтобы
  частица помещалась, а в скриптах ставьте внешний таймаут. Подробнее — раздел 7.6
  [REPORT.md](REPORT.md).
* Память растёт примерно как **0.34 кБ на ячейку решётки для каждого потока** плюс
  около 0.38 кБ на ячейку для матрицы взаимодействия, общей для всех потоков.
  Один поток на решётке 128³ занимает примерно 1.4 ГБ. Точный расчёт и наибольшая
  решётка, влезающая в заданный объём ОЗУ, — в разделе 3 [REPORT.md](REPORT.md).

---

## Состав проекта

| Файл | Содержимое |
|---|---|
| `main.cpp` | цикл по образцам, пул потоков, вывод |
| `options.cpp` | разбор командной строки (argtable3), проверки, интерактивный опрос |
| `dda_func.cpp` | решатель DDA, свёртки, скалярные произведения, поле в дальней зоне |
| `fft_func.cpp` | трёхмерное БПФ |
| `fft_kernel.hpp` | смешанное radix-4/8 ядро 1-D БПФ, используется в `fft_func.cpp` |
| `particle.cpp` | генерация частиц |
| `complex.cpp`, `cvector.cpp` | комплексная арифметика, 3-векторы и матрицы 3×3 |
| `mmuller.cpp` | матрицы Джонса и Мюллера |
| `randomize.cpp` | генерация случайных чисел |
| `function.hpp`, `cplx.h`, `cvector.hpp`, `complex.hpp`, `mmuller.hpp`, `randomize.hpp`, `options.hpp` | объявления |
| `compile` | скрипт сборки |
| `third_party/argtable3/` | вендоренная argtable3 v3.3.1 (BSD-3-Clause), разбор опций |
| `bench/` | отдельный харнесс для 1-D ядра: дельтатесты, точность, A/B-замеры |
| `REPORT.md` | отчёт об оптимизации: профиль, изменения, замеры, отвергнутые идеи |
