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

Version 6.3.

---

## Authors

| | |
|---|---|
| Original program | **Zubko E.S.** — Version 6.3, Copyright (C) 2001–2012 |
| Optimization, parallelization and refactoring | **Shmirko K.A.** — <kshmirko@dvo.ru> |

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

Requires a C++17 compiler: GCC 7+ or clang 5+. No external libraries.

```sh
sh compile          # g++ -O3 -std=c++17 -pthread *.cpp -o dda_fft.out
```

The `compile` script is a plain two-line wrapper; on a system where `g++` is a
different compiler it can be replaced by any equivalent command, e.g.

```sh
g++ -O3 -std=c++17 -pthread *.cpp -o dda_fft.out
```

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
| `fft_kernel.hpp` | mixed-radix 4/8 in-place 1-D FFT kernel used by `fft_func.cpp` |
| `particle.cpp` | particle generation |
| `complex.cpp`, `cvector.cpp` | complex arithmetic, 3-vectors and 3×3 matrices |
| `mmuller.cpp` | Jones and Mueller matrices |
| `randomize.cpp` | random number generation |
| `function.hpp`, `cplx.h`, `cvector.hpp`, `complex.hpp`, `mmuller.hpp`, `randomize.hpp` | declarations |
| `compile` | build script |
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

Версия 6.3.

---

## Авторы

| | |
|---|---|
| Исходная программа | **Зубко Е.С. (Zubko E.S.)** — Version 6.3, Copyright (C) 2001–2012 |
| Оптимизация, распараллеливание и рефакторинг | **Шмирко К.А. (Shmirko K.A.)** — <kshmirko@dvo.ru> |

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

Требуется компилятор с C++17: GCC 7+ или clang 5+. Внешние библиотеки не нужны.

```sh
sh compile          # g++ -O3 -std=c++17 -pthread *.cpp -o dda_fft.out
```

Скрипт `compile` — обёртка из двух строк; на системе, где `g++` означает другой
компилятор, его можно заменить эквивалентной командой, например

```sh
g++ -O3 -std=c++17 -pthread *.cpp -o dda_fft.out
```

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

### Позиционные аргументы

Любой аргумент можно опустить — тогда программа спросит его в интерактивном
режиме, поэтому интерактивный запуск сводится к `./dda_fft.out`. Для скриптов
аргументы позиционные; их смысл:

| № | Аргумент | Примечание |
|---|---|---|
| 1 | `theta_inf` | нижняя граница диапазона фазовых углов, градусы |
| 2 | `theta_up` | верхняя граница диапазона фазовых углов, градусы |
| 3 | `theta_step` | шаг по фазовому углу, градусы |
| 4 | файл `dst` | двоичный файл состояния (текущие средние, см. ниже) |
| 5 | файл `dat` | текстовый файл результата |
| 6 | число образцов | сколько реализаций / ориентаций усреднять |
| 7 | режим ориентации | `0` — фиксированная ориентация, `1` — случайные ориентации (Монте-Карло) |
| 8 | *(режим 0)* `psi` | первый поворот, градусы |
| 9 | *(режим 0)* `iota` | второй поворот, градусы |
| 10 | *(режим 0)* `kappa` | третий поворот, градусы |
| 11 | *(режим 0)* число плоскостей рассеяния | |
| 8 | *(режим 1)* число плоскостей рассеяния | |
| — | тип частицы | см. таблицу ниже; он сдвигает все последующие аргументы |
| — | параметры частицы | зависят от типа (например `m′ m″` для сферы) |
| — | `i_max` `j_max` `l_max` | размеры решётки, степени двойки |
| — | `x` | параметр размера частицы |
| — | модель поляризуемости | `1`…`5`, см. ниже |
| — | `Error_max` | порог невязки, при котором итерации прекращаются |
| — | `Max_Step` | максимальное число итераций |
| — | метод FFT | `0` — Кули–Тьюки (**рекомендуется**), `1` — Радера–Бреннера |
| — | метод решения | `0`…`3`, см. ниже |
| — | `rTau` | только для метода `3` — весовой множитель |

Число фазовых углов равно `(theta_up − theta_inf)/theta_step + 1`.

### Пример: сфера, решётка 64³

```sh
# фазовые углы 0..180 с шагом 10°, фиксированная ориентация, 1 образец,
# одна плоскость рассеяния, сфера, m = 1.5 + 0.01i, 64 x 64 x 64, x = 20,
# Клаузиус–Моссотти, Error_max = 1e-3, 20 шагов, Кули–Тьюки, метод решения 0
./dda_fft.out 0 180 10 out.dst out.dat 1 0 0 0 0 1 1 1.5 0.01 \
              64 64 64 20 1 0.001 20 0 0
```

Для усреднения по 100 случайным ориентациям нужно заменить число образцов и режим
ориентации и указать число плоскостей рассеяния:

```sh
./dda_fft.out 0 180 10 out.dst out.dat 100 1 1 1 1.5 0.01 \
              64 64 64 20 1 0.001 20 0 0
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
DDA_THREADS=4 ./dda_fft.out …
```

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
* Память растёт примерно как **0.34 кБ на ячейку решётки для каждого потока** плюс
  около 0.38 кБ на ячейку для матрицы взаимодействия, общей для всех потоков.
  Один поток на решётке 128³ занимает примерно 1.4 ГБ. Точный расчёт и наибольшая
  решётка, влезающая в заданный объём ОЗУ, — в разделе 3 [REPORT.md](REPORT.md).

---

## Состав проекта

| Файл | Содержимое |
|---|---|
| `main.cpp` | разбор аргументов, цикл по образцам, пул потоков, вывод |
| `dda_func.cpp` | решатель DDA, свёртки, скалярные произведения, поле в дальней зоне |
| `fft_func.cpp` | трёхмерное БПФ |
| `fft_kernel.hpp` | смешанное radix-4/8 ядро 1-D БПФ, используется в `fft_func.cpp` |
| `particle.cpp` | генерация частиц |
| `complex.cpp`, `cvector.cpp` | комплексная арифметика, 3-векторы и матрицы 3×3 |
| `mmuller.cpp` | матрицы Джонса и Мюллера |
| `randomize.cpp` | генерация случайных чисел |
| `function.hpp`, `cplx.h`, `cvector.hpp`, `complex.hpp`, `mmuller.hpp`, `randomize.hpp` | объявления |
| `compile` | скрипт сборки |
| `bench/` | отдельный харнесс для 1-D ядра: дельтатесты, точность, A/B-замеры |
| `REPORT.md` | отчёт об оптимизации: профиль, изменения, замеры, отвергнутые идеи |
