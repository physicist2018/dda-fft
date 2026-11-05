# Оптимизация DDA/FFT кода: отчёт

**Русский** | [English](#english)

Краткое описание программы, авторы и инструкция по сборке и запуску — в
**[README.md](README.md)**. Этот документ посвящён только оптимизации: профилю,
изменениям, замерам и тому, что было проверено и отвергнуто.

Ветка: `perf/stages-0-3` (база — `master`, коммит «Baseline: DDA/FFT source as received»)

Машина, на которой сняты все замеры: **Apple M1, 8 ядер (4P + 4E), 8 ГБ ОЗУ, Apple clang 16**.
Компиляция: `sh compile`, то есть `cc -O2 -std=c11 -c third_party/argtable3/argtable3.c` и
`g++ -O3 -std=c++17 -pthread *.cpp third_party/argtable3/argtable3.o -o dda_fft.out`.

Версия программы — **7.2**; изменение и оптимизация: **К.А. Шмирко** <kshmirko@dvo.ru>, 2026.

---

## 1. Краткое резюме

| Показатель | Было | Стало | Ускорение |
|---|---|---|---|
| 64³, 19 углов, 20 шагов, 1 образец | 29.6 с | **8.85 с** | **3.34×** |
| 128³, 19 углов, 3 шага, метод 2 | 120.3 с | **35.7 с** | **3.37×** |
| 32³, 19 углов, 40 шагов | 6.33 с | **1.83 с** | 3.46× |
| 64³, 181 угол, 3 шага | 12.2 с | **4.03 с** | 3.02× |
| **64³, 20 реализаций × 20 шагов (`W=7`)** | **576.4 с** | **76.6 с** | **7.5×** |

Три источника ускорения (этапы 0–4):

1. **Бит-инверсия в FFT** — 2.85× (главный вклад).
2. **Непрерывные блоки памяти + пул рабочих буферов + паддинг строк** — ещё ~1.17×.
3. **Параллельные вычисления по реализациям** — 2.2× (на M1; на машинах с большей полосой памяти и большим числом P-ядер должно быть лучше).

Плюс этап 6 — собственное radix-4/8 ядро 1-D FFT (раздел 4): на той же задаче 64³,
20 шагов, один поток — **5.67 с → 4.00 с** (ещё 1.42×), а относительно исходного кода
**18.98 с → 4.00 с**, то есть **4.75× суммарно**. Замеры этапа 6 сняты в отдельной сессии
(см. оговорку в начале раздела 5): сравнивать между разделами надо отношения, а не
абсолютные секунды.

Память выросла на **1 %**. Через этапы 0–4 результаты на детерминированных задачах были
**побитово идентичны** исходным; после этапа 6 побитовость исчезает — другой порядок
суммирования в `float` даёт расхождение порядка 10⁻⁴ на самых слабых значимых элементах
матрицы (раздел 4, «Этап 6»).

Плюс исправлен **дефект чтения неинициализированной памяти** в итерационном процессе.

---

## 2. Исходное состояние

Проект — реализация Discrete Dipole Approximation: рассеяние света на ансамбле диполей.
Ядро — итерационное решение системы `(I − αA)P = αE`, где произведение на `A` считается
трёхмерной свёрткой через FFT. Размерность задачи задаётся решёткой `i_max × j_max × l_max`,
реально используются кубы 64³ и 128³; внутри свёртка работает на массиве, дополненном
нулями в 2 раза по каждой оси (то есть 8·N ячеек).

Структура вызовов:

```
main.cpp          цикл реализаций / ориентаций
  particle.cpp    генерация формы частицы
  dda_func.cpp    DDA_FFT: итерационное решение
     Convolution  -> FFT_3D / iFFT_3D  (fft_func.cpp)
     LinearComb, Copying, ldcScalMul, ldError   (покомпонентные циклы)
     CabsCext, FarField
```

---

## 3. Где на самом деле уходит время

Профиль снят инструментированием текущей сборки (64³, 20 шагов, 1 образец, `W=1`, 9.08 с):

| Компонент | Время | Доля |
|---|---|---|
| **FFT** (`FFT_3D` + `iFFT_3D`) | 7.27 с | **80 %** |
| Точечное умножение в свёртке | 0.61 с | 7 % |
| `CopyingBack` | 0.31 с | 3 % |
| Остальное (элементные циклы, `FarField`, шаг 0) | ~0.9 с | 10 % |

Внутри FFT: **1-D ядро ≈ 63 %**, сборка/разборка строк и умножение на твидал ≈ 37 %.

Память: **673 байта на ячейку** решётки всего, из них

| Массив | Байт/ячейка | Комментарий |
|---|---|---|
| `cmA` (ДПФ ядра взаимодействия) | 384 | общая для всех потоков |
| `cvWork` (дополненный буфер свёртки) | 192 | на поток |
| `cvP`, `cvP_min`, `cvR`, `cvS` | 96 | на поток |
| `cvG`, `cvH` (только метод 2) | 48 | на поток |
| `cAlpha` | 1 | |

Отсюда: 128³ ≈ 1.4 ГБ при одном потоке; 256³ (11.3 ГБ) в 8 ГБ не влезает в принципе.

---

## 4. Что сделано

### Этап 0 — сборка (коммит `df8a07a`)

* `compile` переведён на `-std=c++17 -pthread`. Раньше Apple clang собирал файлы как **C++98**;
  без C++17 недоступны `std::thread` и лямбды.
* Убраны 12 мест с сужением `double → float` в списках инициализации (в `dda_func.cpp` и
  `particle.cpp`) — под C++11+ это ошибка компиляции, причём **и в GCC тоже**, так что правка
  нужна и для Intel/AMD.
* Исправлены 4 предупреждения `printf("%d", time_t/clock_t)` → `%ld` (переносимость).
* `-march=native` сознательно **не** используется: замерено, что он даёт 0 %, а бинарник
  остаётся переносимым.

### Этап 1 — бит-инверсия (`a3cf165`, `fft_func.cpp`) — **2.85×**

```c
/* было: malloc + free + pow() на КАЖДЫЙ элемент */
bits = malloc(number_of_bits*sizeof(short int));
...
(*number_new) += (int)pow(2., (double)i) * bits[i];

/* стало */
for (i = 0; i < number_of_bits; i++)
    r = (r << 1) | ((number_old >> i) & 1);
```

Функция вызывается один раз на точку FFT, то есть для преобразования 128³ — **6.3 млн раз**
на одно трёхмерное преобразование, и дважды за итерацию (прямое + обратное).
Микробенчмарк: 69 нс/вызов → 3.5 нс/вызов. На реальной задаче: 29.6 с → 10.4 с, вывод
**побитово идентичен**.

### Этап 2 — память и пул буферов (`bc67c5c`) — ещё 1.17×

* `alloc_field3()`: массив `ni×nj×nl` выделяется **одним блоком** (3 `malloc` вместо ~33 000
  на вызов), указатели-строки нарезаются внутрь блока. API `cvector***` сохранён, поэтому
  все существующие циклы не тронуты.
* `dda_workspace` — семь рабочих массивов живут весь прогон, а не выделяются/освобождаются
  дважды на каждую реализацию. Это одновременно устранило повторные page-fault-штормы
  и стало предпосылкой для этапа 4.

**Важная находка.** Просто непрерывная раскладка дала **регресс 6 %**: размеры решётки —
степени двойки, поэтому все строки попадали на кратные степени двойки смещения, и три
прохода FFT «долбили» одни и те же множества кэша. Добавлен `FIELD_PAD = 1` (паддинг строки
в один элемент, 24 байта) — стало **+17 % вместо −6 %**, ценой +0.8 % памяти. Проверялись
PAD = 1, 2, 4, 8 — единицы достаточно.

### Исправление дефекта (`cc64d5c`)

В конце `DDA_FFT` есть откат к «лучшей» итерации:

```c
if (Error > Error_min) cvP = cvP_min;
```

но `cvP_min` заполнялся **только если какая-то итерация улучшила невязку**. Если процесс ни
разу не улучшился (типично для метода простой итерации с большим весовым множителем, а также
возможно для остальных методов на трудных задачах), откат копировал в дипольные моменты
**неинициализированную память**. Проявление: один и тот же бинарник с одними и теми же
аргументами давал `Cabs` = 4.5e7, 2.0e6 и 5.2e10 в трёх прогонах подряд.

Исправление — итерация шага 0 сохраняется сразу после вычисления. Поведение не меняется,
если хоть одна итерация улучшает невязку (там `cvP_min` всё равно перезаписывается).

Следствие для сверки со старыми результатами: **выход метода 3 больше не совпадает с тем,
что давал исходный код, и не должен** — раньше там был мусор из кучи (разный от прогона к
прогону), теперь это итерация шага 0. Для проверки: при расходящейся итерации результат
`Max_Step=5` теперь побитово равен результату `Max_Step=0`, а повторные прогоны совпадают
между собой. Методы 0, 1 и 2 от исправления не меняются (проверено побитово на 32³, 64³,
128³ и на неквадратной решётке 32×64×32).

### Этап 4 — параллельные вычисления (`5f81e60`, `ea7530d`, `6b6f799`) — 2.2×

Реализации независимы, поэтому цикл по ним распараллелен:

* **4a**: `rand()` хранит состояние в глобальной переменной и не потокобезопасен. Добавлен
  потоковый ГПСЧ (`thread_local` xorshift64*) за `dda_rand()`, возвращающий тот же 31-битный
  диапазон; все вызовы в `particle.cpp`, `dda_func.cpp`, `main.cpp` переведены на него.
* **4b**: в `dda_workspace` добавлен флаг `iVerbose` — при нескольких потоках пошаговые
  строки невязки не печатаются (иначе они перемешиваются), вместо них одна сводка на
  готовую реализацию.
* **4c**: тело цикла реализаций вынесено в лямбду; внутри неё локальные имена привязаны к
  приватным массивам воркера, поэтому **сам код цикла остался прежним** — изменились только
  аккумуляторы и печать. Воркеры берут индексы из `std::atomic`; накопление сумм и запись
  `.dst` — под мьютексом один раз на реализацию (файл сохранил смысл, прогон можно прервать
  и продолжить).

Приватно у воркера: карта занятости, поляризуемости, геометрия рассеяния, поля `Einc`/`Esc`,
`dda_workspace`. Общие только на чтение: `cmA`, таблицы FFT, `ksc_ini`/`vNormal_ini`.

### Этап 6 — собственное radix-4/8 ядро 1-D FFT (ветка `perf/fft-radix48`)

Это последний пункт из раздела 11 предыдущей версии отчёта. Старое ядро (`fft_func.cpp`) —
radix-2 DIT, которое (а) пересчитывало бит-инверсию на каждый вызов и (б) после каждой
стадии, начиная со второй, копировало всю строку из одного массива в другой. Для строки из
N точек это `log2(N)−1` полных копий на преобразование плюс отдельный проход перестановки,
и именно трафик памяти, а не арифметика, ограничивал ядро.

Что сделано (`fft_kernel.hpp`, header-only; старые `FFT`/`iFFT` оставлены в `fft_func.cpp`
как резервный путь):

* План на каждый размер N строится один раз и дальше только читается: смешанная
  digit-reversal перестановка под выбранную последовательность радиксов и полная таблица
  `W_N^k` из N элементов, посчитанная напрямую через `cos`/`sin` в `double` и приведённая
  к `float`. Планы кешируются (до 16 размеров) и строятся под мьютексом, потому что
  `FFT_3D` вызывается из нескольких воркеров одновременно.
* Стадии идут **на месте**, в выходном буфере: нет пинг-понга между стадиями и нет
  промежуточного `cvWork2`. Вход не модифицируется, но `in` и `out` не могут совпадать:
  строка сразу загружается в переставленном порядке (`out[k]=in[perm[k]]`). Сам драйвер
  по-прежнему собирает строку в `cvWork` и разбирает результат из `cvResult`.
* Radix-8 и radix-4 бабочки (radix-8 собирается из двух radix-4); три компоненты `x, y, z`
  обрабатываются вместе. При `i=0` умножений на твиды нет; внутри radix-8 твиды `±i`
  бесплатны, а `±(1±i)/√2` — по два умножения на константу 1/√2 на компоненту; на
  остальных элементах — обычные комплексные умножения на табличные твиды.
* Разложение `log2(N)` выбирается по таблице (`8·8`, `8·4`, `4·4`, `8·4·4`, …), крупный
  radix — в первые, самые локальные стадии; radix-2 остаётся только для N=2.
* Ядро подключается только для степеней двойки: при `N = 2^k` драйвер берёт план, иначе
  (или при сборке с `-DFFT_KERNEL_LEGACY`) работает прежнее radix-2 ядро. Обратное
  преобразование считается тем же ядром по схеме `conj → DFT → conj`, без множителя `1/N`,
  как и старый `iFFT`.

Контракт сохранён: ядро по-прежнему считает **`DFT_N` в естественном порядке**. «Половинная»
конвенция живёт не в ядре, а в `FFT_3D`: дополненная строка длины `2N` имеет только `N`
ненулевых отсчётов, поэтому `X[2l] = DFT_N(x)[l]`, `X[2l+1] = DFT_N(x·W_{2N}^j)[l]`, и драйвер
перемежает чётные/нечётные бины. Это подтверждено дельтатестом: старое ядро — обычное
`DFT_N`, а составное преобразование строки — полное `DFT_2N` дополненной линии.

Замеры (та же машина, что и в остальных разделах):

| Что | Было | Стало | Ускорение |
|---|---|---|---|
| 1-D ядро, N=16 | 131 нс/строка | 77 | 1.71× |
| 1-D ядро, N=32 | 282 | 184 | 1.53× |
| 1-D ядро, N=64 | 632 | 390 | 1.62× |
| 1-D ядро, N=128 | 1459 | 940 | 1.55× |
| 1-D ядро, N=256 | 3299 | 2121 | 1.56× |
| 64³, 20 шагов, 1 образец, `W=1` | 5.67 с | **4.00 с** | **1.42×** |
| 128³, 5 шагов, 1 образец, `W=1` | 16.75 с | **11.56 с** | **1.45×** |
| 64³, 4 образца × 5 шагов, `W=4` | 2.23 с | **1.58 с** | **1.41×** |
| 64³, исходный код → эта ветка | 18.98 с | **4.00 с** | **4.75×** |

Абсолютные наносекунды и отношения в таблице — из сессии этапа 6. Повторный прогон
`compare` на той же машине даёт 1.5–1.6× для N=16…256, абсолютные значения плавают на
несколько процентов от загрузки машины.

Оценка из прежнего раздела 11 (10–15 % сквозного времени) оказалась заниженной: ядро
ускорилось не в 1.2–1.5, а в 1.5–1.7 раза, что при доле FFT ≈80 % даёт 1.42×. Это ровно
предел Амдала для ускорения ядра 1.6×: `0.8/1.6 + 0.2 = 0.70`.

**Точность.** Само ядро считает точнее старого (ошибка относительно `double`-ДПФ на N=128:
1.1·10⁻⁶ → 8.9·10⁻⁸), потому что таблица твидов считается напрямую, а не рекуррентным
умножением. Но результат перестал быть **побитово** идентичным: `real` здесь `float`, и
другой порядок суммирования даёт другой шум округления. Проверено:

* дельтатест нового ядра по всем позициям, N=4…256: ошибка ≤ 7.3·10⁻⁸;
* `float` master против `double` master и `float` new против `double` master на 64³:
  худшие значимые элементы расходятся на 1.2·10⁻⁴ и 1.9·10⁻⁴ — то есть оба `float`-бинарника
  одинаково далеки от точного ответа, а между собой — 8.8·10⁻⁵; расхождение растёт с числом
  итераций (1.2·10⁻⁶ за 1 шаг → 5·10⁻⁵ за 10), как и положено накоплению округления;
* сборка с `-DFFT_KERNEL_LEGACY` (старое ядро) воспроизводит master **побитово** — это
  регрессионный тест на то, что переделка драйвера ничего не изменила.

Проверки и замеры воспроизводятся харнессом `bench/` (`sh bench/build.sh`; режимы `delta`,
`newdelta`, `accuracy`, `compare`, `half`, `speed`, `roundtrip`), а сравнение выходных `.dat`
между сборками — скриптом `bench/cmp_dat.py`. В основной бинарник харнесс не входит:
`compile` собирает `*.cpp` из корня, а харнесс лежит в подкаталоге.

---

## 5. Результаты

*Оговорка.* Таблицы этого раздела — замеры этапов 0–4, снятые в одной сессии. Замеры
этапа 6 (раздел 4) сняты позже, в другой сессии и на более свободной машине, поэтому
абсолютные секунды между разделами не сравнимы — сравнимы только отношения. Так, в этой
сессии `master` даёт на 64³/20 шагов 5.67 с (а не 8.85 с, как в таблице ниже), и новый
этап 6 даёт 4.00 с, то есть те же 1.42×.

### Сквозные замеры

| Конфигурация | Было | Стало | Ускорение |
|---|---|---|---|
| 32³, 19 углов, 40 шагов | 6.33 с | 1.83 с | 3.46× |
| 64³, 19 углов, 20 шагов | 29.6 с | 8.85 с | 3.34× |
| 64³, 181 угол, 3 шага | 12.2 с | 4.03 с | 3.02× |
| 128³, 19 углов, 3 шага, метод 2 | 120.3 с | 35.7 с | 3.37× |
| 64³, 20 реализаций × 20 шагов, `W=1` | 576.4 с | 171.1 с | 3.37× |
| 64³, 20 реализаций × 20 шагов, **`W=7`** | 576.4 с | **76.6 с** | **7.5×** |

### Масштабирование по потокам (64³, 20 реализаций × 20 шагов)

| W | 1 | 2 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| время, с | 171.1 | 120.9 | 100.3 | 88.7 | 84.6 | **76.6** | 80.4 |
| ускорение | 1.00 | 1.42 | 1.71 | 1.93 | 2.02 | **2.23** | 2.13 |

128³ (4 реализации, 2 шага): `W=1` 105.2 с, `W=2` 74.5 с (**1.41×**), `W=4` 78.7 с.

**Почему на 128³ оптимум — 2 потока.** Рабочее множество воркера там 674 МБ; два воркера
уже насыщают полосу памяти M1 (~68 ГБ/с), третий и четвёртый не помогают. На 64³ множество
втрое меньше, поэтому выгодно 7 потоков (8-й конкурирует с системой).

**Накладных расходов от потоков нет:** 4 потока проходят ту же работу ровно за столько же,
сколько 4 независимых однопоточных процесса (100.3 с против 100.5 с).

### Память

Пиковый RSS на 128³: 1.508 → **1.524 ГБ** (+1 %). Рост — только паддинг строк.

---

## 6. Что проверено и отвергнуто (с доказательствами)

Этот раздел важен не меньше сделанного: он фиксирует, куда не надо возвращаться.

### 6.1. Компактизация полей по активным ячейкам — только ёмкость, не скорость

Сфера занимает 52 % бокса, поэтому хранение `cvP/cvR/cvS/cvP_min` только по активным ячейкам
дало бы −16 % памяти на воркера. Но на 64³/128³ память **не является ограничением** (1.4 ГБ
при 8 ГБ), а сетку 256³ это всё равно не открывает: одна `cmA` для неё — 6.4 ГБ. Выигрыша по
скорости нет: элементные циклы обслуживаются в тех ~10 % времени, которые остались после FFT.

### 6.2. Уменьшить `cmA` нельзя — проверено численно

Гипотеза: ядро имеет вид `A(r) = f(|r|)·I + g(|r|)·dyad(r)`, поэтому его фурье-образ должен
быть `Â(k) = A(k)·I + B(k)·k̂k̂` — всего 2 комплексных числа вместо 6 (48 → 16 байт на ячейку,
−36 % памяти и −67 % трафика свёртки). Написан диагностический расчёт на реальном
`CreateAMatrix` + `FFT_3D`:

```
grid 16³: max relative deviation from A*I+B*kk^ = 2.19e+01
          (только k внутри |k|d0 < π/2)         = 8.67e-01
grid 32³: max relative deviation                = 1.38e+01
          (только k внутри |k|d0 < π/2)         = 1.04e+00
```

Отклонение ~100 % даже в центре зоны Бриллюэна: кубическая анизотропия дискретизированного
на решётке ядра сравнима с изотропной частью. Тождество **не выполняется**, идея закрыта.

### 6.3. Отказ от зануления неиспользуемой части `cvWork` — 0 %

Для режима 0 FFT действительно не читает 7/8 дополненного буфера (проверено по границам
циклов всех трёх проходов). Убрал зануление: **8.85 с против 8.81 с**, результат побитово
идентичен. Цикл зануления и так дёшев и дружелюбен к кэшу.

### 6.4. Этап 5, пункты 1–2: SoA и пакетная обработка строк — 0.37×…1.01×

Собрано шесть вариантов 1-D ядра, каждый замерен отдельно:

| Вариант | нс/строка | Ускорение |
|---|---|---|
| A — текущий код (AoS `cvector`) | 1209 | 1.00× |
| B — SoA по плоскостям (6 массивов) | 2581 | **0.37×** |
| C — AoS, 6 подряд идущих `float` в бабочке | 1220 | 1.00× |
| E — замена целочисленного деления `i*n/iStep` на умножение | 1193 | 1.01× |
| всё то же с `-march=native` | — | 1.00× |
| D — батч из 4 строк, раскладка `[l][4]`, вектор по строкам | 1132 | **0.90×** |
| D + принудительная векторизация (`#pragma clang loop vectorize`) | 1130 | 0.91× |

Причины: при строке в 64 элемента разнесение по плоскостям ломает локальность (бабочка
трогает 6 участков вместо двух рядом лежащих блоков); в большинстве стадий внутренний цикл
имеет 1–2 итерации, а параллелизм лежит по индексу строки, где данные идут с шагом —
компилятор прямо сообщает `cost-model indicates that vectorization is not beneficial`.

### 6.5. pocketfft — неприменим (из-за структуры задачи, не из-за качества)

Библиотека взята из `github.com/mreineck/pocketfft` (BSD-3, header-only, портируема).

*Замечание о целостности файла*: строки вида `cb.i=twai*t4.r twbi*t3.r;` выглядят как
повреждение, но это идиома автора — знак приходит из аргумента макроса
(`POCKETFFT_PARTSTEP5a(...,+tw1i,+tw2i)`). Файл корректен, собирается и верно считает для
всех размеров (проверено n = 2, 4, 5, 7, 8, 10, 11).

Замеры:

| Способ интеграции | 128³, 3 компоненты |
|---|---|
| наш `FFT_3D` (как есть) | **71 мс** |
| pocketfft `c2c` по всему дополненному объёму | 102 мс |

pocketfft **в 1.4 раза медленнее**, потому что наш код знает про нулевой паддинг: на первом
проходе он читает только значащий блок (64³ из 128³) и разворачивает результат, то есть
делает примерно **вчетверо меньше работы**. Библиотека общего назначения вынуждена
преобразовывать весь объём.

Интеграция только внутреннего ядра тоже невозможна: дельта-тест показывает, что «1-D FFT»
в этом коде — это N-точечное преобразование, а полный спектр дополненной строки длины 2N
собирается из двух таких преобразований (отсюда и перемежение `2*l` / `2*l+1` у
вызывающего кода):

```
дельта в позиции 1:
  текущее ядро:   +0.92-0.38i  +0.71-0.71i  +0.38-0.92i     шаг π/8  → 2n-точечное ДПФ
  pocketfft:      +0.71-0.71i  +0.00-1.00i  -0.71-0.71i     шаг π/4  → n-точечное ДПФ
```

Корректная замена потребовала бы считать 2n-точечное преобразование — вдвое больше работы,
что съедает весь выигрыш ядра (pocketfft быстрее «за единицу работы» лишь в 1.17× при n=64
и 1.48× при n=128).

**Осторожно:** в ходе этой проверки прототип с подменённым ядром показал «+3 % на 64³ и
+19 % на 128³», но эти числа **недействительны** — он считал другое преобразование. Именно
дельта-тест это и выявил. Не использовать эти цифры.

Заголовок библиотеки (`pocketfft_hdronly.h`, BSD-3, 121769 байт, sha256
`518a95c5…44c86d`) в репозиторий не входит и с проектом не распространяется — он был нужен
только для этой проверки. Если понадобится повторить замеры, его можно взять из
`github.com/mreineck/pocketfft`, ветка `cpp`.

---

## 7. Найденные дефекты

### 7.1. Чтение неинициализированной памяти в `cvP_min` — **исправлено** (`cc64d5c`)

См. раздел 4. Проявление: недетерминированный результат метода простой итерации
(`Cabs` = 4.5e7 / 2.0e6 / 5.2e10 в трёх прогонах).

### 7.2. Режим `iMode_FFT = 1` (Rader–Brenner) асимптотически плох — **не исправлено**

Внутри 1-D `FFT` для этого режима стоит двойной цикл `for(i < N/2) for(j < N/2)`, то есть
**O(L²) на строку** и O(N²) на трёхмерное преобразование. Уже при 32³ он медленнее режима 0
(1.36 с против 1.00 с), и разрыв растёт с размером.

**Рекомендация: не использовать `iMode_FFT = 1`.** Режим 0 (Cooley–Tukey) — единственный
практичный. Код режима 1 оставлен без изменений, чтобы не менять поведение.

### 7.3. Сборка по умолчанию как C++98 — **исправлено** (`df8a07a`)

### 7.4. `DDA_FFT` ранее не был реентерабельным — следствие исправлено попутно

До этапа 2 рабочие массивы выделялись/освобождались внутри `DDA_FFT`, что само по себе
исключало параллелизм. Теперь состояние вынесено в `dda_workspace`.

### 7.5. `FourierPhase(1)` пишет за границу массива — **не исправлено**

Найдено при отладке харнесса этапа 6. `FourierPhase(i_max, cW)` безусловно заполняет
`cW[0]` и `cW[1]`, а таблицы `cW_2_*` выделяются размером `i_max/2`. При размерности
решётки 2 это `FourierPhase(1)` в массив из одного элемента — запись на 8 байт за границу
кучи. Достижимо только для решётки с размером 2 (физически бессмысленной), и только в
старом коде; новое ядро эту таблицу не использует, но `main.cpp` её всё равно строит.
Оставлено без правки, чтобы не менять поведение за рамками задачи.

### 7.6. Зависание генерации частиц — **не исправлено**

Найдено при регрессе нового интерфейса: все 14 типов целей сверялись с эталонной сборкой
побайтово, и на двух наборах параметров **зависают обе сборки** — и старая (позиционные
аргументы), и новая. То есть дефект существовал и раньше и с разбором командной строки
не связан.

**`iGenerateHedgehog` (`particle.cpp`)** — та же схема размещения, что в
`iGenerateAnthrax` и `iGenerateDeadAnthrax`:

```c
while(m<iWormsNumber)
{   cl0.r=-1.;
    while(cl0.r<0.)      /* цилиндр не влез в матрицу — тянем заново */
    {   ...
        cl0=worm(...,length_worm,radius_worm,R2);
    }
    for(j=0;j<m;j++) { l=CrossTwoWorms(cl0,worms[j]); if(l==1) break; }
    if(l==0) { worms[m]=cl0; m++; }
}
```

Оба цикла бесконечны, если допустимого положения нет: внутренний — когда мономер
(`--monomer-radius` вместе с `--monomer-length-ratio`) не помещается в матрицу, внешний —
когда `iWormsNumber`, вычисленное из `--packing`, больше, чем можно разместить без
пересечений. Проверки числа попыток нет, поэтому программа просто не возвращается из
генерации.

**`iGenerateHierarchicalEllipsoid`** — поиск затравочной ячейки:

```c
do
{   iFlag=0;
    do
    {   v0 = случайная ячейка матрицы;
        i0 = cAlpha[(int)v0.x][(int)v0.y][(int)v0.z];
    }
    while(i0!=1);        /* нужна ячейка именно первого материала */
    ...  /* отбраковка уже выбранных seed'ов */
}
while(iFlag==1);
vSeeds[m]=v0;
```

Внутренний цикл не завершится, если в текущем уровне **нет ни одной ячейки со значением 1**
(или их меньше, чем `iNoSeeds`). Типично для тесной матрицы: например, у эллипсоида с малыми
полуосями на решётке 8³ внутренняя область оказывается пустой.

**Воспроизведение** (зависает и старая сборка с теми же параметрами):

```sh
# тип 7: слишком плотная упаковка для 8³
./dda_fft.out --target hedgehog --packing 0.3 --monomer-radius 1 \
    --monomer-length-ratio 1 --material 1.5 0.01 --grid 8 --size-parameter 2

# тип 13: у эллипсоида 0.8 x 0.6 на 8³ не остаётся ячеек первого материала
./dda_fft.out --target hierarchical --aspect 0.8 0.6 --levels 2 \
    --material 1.5 0.01 --material 1.33 0.0 --seeds 111 222 \
    --grid 8 --size-parameter 2
```

С параметрами, при которых места хватает, оба типа считаются нормально: тип 7 —
`--packing 0.05 --monomer-radius 0.5 --grid 16`, тип 13 — `--aspect 1 1 --levels 2 --grid 16`
(на этих наборах результаты новых и старых сборок совпали побайтово).

**Почему не исправлено.** Правильная правка — ограничить число попыток и вернуть ошибку или
частично заполненную частицу. Это меняет геометрию и, следовательно, числа, то есть требует
отдельной задачи с проверкой физики. До исправления: подбирать `--packing`/`--monomer-*` под
объём матрицы, для типа 13 брать достаточно крупную решётку, а в скриптах ставить внешний
таймаут.

---

## 8. Как пользоваться

### Сборка

```sh
sh compile
```

Требуется компилятор с C++17: GCC 7+ либо clang 5+. На Linux/Intel/AMD подойдёт системный
`g++`; отдельные пакеты не нужны: argtable3 вендорена в `third_party/argtable3/`
(BSD-3-Clause) и собирается тем же `compile` как C.

### Число потоков

```sh
./dda_fft.out --threads 4 <опции>              # явно
./dda_fft.out <опции>                          # по умолчанию: по числу ядер
```

Переменная окружения `DDA_THREADS` учитывается, если `--threads` не задан.

По умолчанию берётся число ядер, но не больше числа реализаций и не больше, чем позволяет
память (70 % физической ОЗУ минус общая `A`-матрица). О выборе печатается строка

```
 Worker threads: N
```

При `N = 1` вывод полностью прежний (пошаговые невязки); при `N > 1` — одна сводка на
готовую реализацию.

**Подбор числа потоков зависит от размера решётки**, потому что ограничение — полоса памяти,
а не ядра. На M1/8 ГБ измеренный оптимум: 64³ → 7 потоков, 128³ → 2 потока. На других
машинах оптимум другой, его стоит подобрать одним прогоном с 1, 2, 4, 8 потоками.

### Запуск

```sh
# 64³ сфера, m = 1.5 + 0.01i, 19 углов, 20 шагов
./dda_fft.out --threads 7 --target sphere --material 1.5 0.01 --grid 64 \
              --size-parameter 20
```

Метод FFT и метод решения задаются опциями `--fft` и `--solver`.

---

## 9. Переносимость (Intel / AMD / Apple Silicon)

* Код правок архитектурно-независим: ни интринсиков, ни `-march=native`.
  Синтаксис под x86-64 проверен (`clang++ --target=x86_64-apple-darwin -fsyntax-only`).
* `FIELD_PAD = 1` лечит конфликты по множествам кэша и на x86 (4 КБ страницы,
  64-байтовые линии): шаг строки перестаёт быть кратным степени двойки.
* **`long double` — не проблема для x86.** На Apple Silicon он равен `double` (8 байт), на
  x86 — 80 бит (x87, без векторизации). Замерена доля трёх функций-сумматоров
  (`ldcScalMul`, `ldScalMulConj`, `ldError`): **0.075 с из 9.5 с = 0.8 %**. Даже при
  3–4-кратном замедлении x87 это +2–3 % на Intel. Трогать не нужно.
* **Результаты не будут побитово совпадать между arm64 и x86** — из-за разной ширины
  `long double` в аккумуляторах. Это свойство кода было и раньше; в пределах одной
  платформы всё детерминировано.

---

## 10. Воспроизведение замеров

Общий вид команды:

```sh
./dda_fft.out --state out.dst --output out.dat --samples <N> \
    --orientation <fixed|mc> --target sphere --material <mre> <mi> \
    --grid <i> <j> <l> --size-parameter <x> --polarizability <1..5> \
    --tolerance <Error_max> --max-steps <Max_Step> --fft <0|1> --solver <0..3>
```

Замеры из отчёта (19 углов; сфера; `m = 1.5+0.01i`; `x = 20`; поляризуемость 1):

```sh
# 64³, 20 шагов, метод 0
./dda_fft.out --state o.dst --output o.dat --target sphere --material 1.5 0.01 \
    --grid 64 64 64 --size-parameter 20 --tolerance 1e-3 --max-steps 20 --fft 0 --solver 0

# 128³, 3 шага, метод 2
./dda_fft.out --state o.dst --output o.dat --target sphere --material 1.5 0.01 \
    --grid 128 --size-parameter 40 --tolerance 1e-3 --max-steps 3 --fft 0 --solver 2

# 64³, 20 реализаций × 20 шагов (главный «пользовательский» сценарий)
./dda_fft.out --threads 7 --state o.dst --output o.dat --samples 20 \
    --target sphere --material 1.5 0.01 --grid 64 --size-parameter 20 \
    --tolerance 1e-3 --max-steps 20 --fft 0 --solver 0
```

Проверки корректности, использованные в работе:

```sh
# побитовая идентичность результата
cmp new.dat old.dat

# потокобезопасность
cc -O1 -std=c11 -c third_party/argtable3/argtable3.c -o argtable3.o
g++ -O1 -g -std=c++17 -pthread -fsanitize=thread  -o t *.cpp argtable3.o && ./t --threads 4 ...
g++ -O1 -g -std=c++17 -pthread -fsanitize=address -o a *.cpp argtable3.o && ./a --threads 4 ...
```

Проверено: методы 0/1/2 (и 3 — см. ниже) и оба режима FFT; решётки 32³, 64³, 128³ и
неквадратная 32×64×32; фиксированные и случайные ориентации; ThreadSanitizer и
AddressSanitizer — чисто.

---

## 11. Что осталось и почему не сделано

Пункт про собственное radix-4/8 ядро — единственный оставшийся рычаг — закрыт этапом 6.
После него FFT-свёртка упирается в предел Амдала: при доле FFT ≈80 % и ускорении ядра 1.6×
остальные 20 % (точечное умножение, сборка/разборка строк, `CopyingBack`, шаг 0)
ограничивают выигрыш величиной ≈1.43×, и он уже достигнут.

Мелкие резервы — совместить сборку строки с перестановкой, вынести рабочие буферы в
`dda_workspace` вместо `malloc` на каждый трёхмерный вызов — дают единицы процентов и не
окупают усложнения.

Дальнейший рост возможен только сменой алгоритма (например, многомасштабные/иерархические
методы вместо FFT-свёртки на больших решётках), а не оптимизацией текущей схемы.

---

## 12. История изменений

```
af63403  Baseline: DDA/FFT source as received
df8a07a  Stage 0: build with C++17, fix narrowing inits and printf format specifiers
a3cf165  Stage 1: replace malloc/pow bit-reversal with shift-based version
bc67c5c  Stage 2: contiguous field storage and a reusable solver workspace
cc64d5c  Fix: keep the step-0 iterate in cvP_min
5f81e60  Stage 4a: per-thread random number generator
ea7530d  Stage 4b: optional quiet mode for the solver
6b6f799  Stage 4c: solve the samples in parallel
(perf/fft-radix48)  Stage 6: mixed-radix 4/8 in-place 1-D FFT kernel
4eae0b7  Vendor argtable3 3.3.1 for command-line parsing
a0f77ac  Replace positional arguments with named command-line options
```

Суммарно относительно базы: **8 файлов, +688 / −388 строк** для этапов 0–4; этап 6 добавляет
`fft_kernel.hpp` и правку `fft_func.cpp` (+67 / −46). Изменения не трогают физику расчёта —
только сборку, раскладку памяти, способ вычисления FFT и порядок обхода независимых
реализаций.

Баннер программы в `main.cpp` обновлён до версии **7.1** с указанием ответственного за
изменение и оптимизацию (`K.A.Shmirko`, <kshmirko@dvo.ru>, 2026).

Ветка `feature/cli-argtable3` (тег **v7.2**) переводит настройку программы на именованные
опции командной строки: разбор вынесен в `options.cpp`/`options.hpp`, позиционные аргументы
удалены, добавлены проверки диапазонов, `--seed` и `--restart`. Версия в баннере, в
`DDA_VERSION` и в документации — **7.2**.

Из стороннего кода в репозитории только вендоренная argtable3 v3.3.1
(`third_party/argtable3/`, BSD-3-Clause); заголовок pocketfft из раздела 6.5 в поставку
не входит.

---

<a id="english"></a>

# Optimising the DDA/FFT code: report

A short description of the program, the authors, and build/run instructions are
in **[README.md](README.md)**. This document covers only the optimisation work:
the profile, the changes, the measurements, and what was tried and rejected.

Branch: `perf/stages-0-3` (base — `master`, commit "Baseline: DDA/FFT source as received")

All measurements were taken on: **Apple M1, 8 cores (4P + 4E), 8 GB RAM, Apple clang 16**.
Build: `sh compile`, i.e. `cc -O2 -std=c11 -c third_party/argtable3/argtable3.c` and
`g++ -O3 -std=c++17 -pthread *.cpp third_party/argtable3/argtable3.o -o dda_fft.out`.

Program version — **7.2**; modification and optimisation: **K.A. Shmirko** <kshmirko@dvo.ru>, 2026.

---

## 1. Summary

| Case | Before | After | Speed-up |
|---|---|---|---|
| 64³, 19 angles, 20 steps, 1 sample | 29.6 s | **8.85 s** | **3.34×** |
| 128³, 19 angles, 3 steps, solver 2 | 120.3 s | **35.7 s** | **3.37×** |
| 32³, 19 angles, 40 steps | 6.33 s | **1.83 s** | 3.46× |
| 64³, 181 angles, 3 steps | 12.2 s | **4.03 s** | 3.02× |
| **64³, 20 samples × 20 steps (`W=7`)** | **576.4 s** | **76.6 s** | **7.5×** |

Three sources of speed-up (stages 0–4):

1. **Bit reversal in the FFT** — 2.85× (the main contribution).
2. **Contiguous storage + a pool of workspaces + row padding** — another ~1.17×.
3. **Solving the samples in parallel** — 2.2× (on an M1; machines with more memory
   bandwidth and more P-cores should do better).

On top of that, stage 6 replaced the 1-D kernel with an own radix-4/8 transform
(section 4): on the same 64³, 20-step, single-thread case **5.67 s → 4.00 s** (a further
1.42×), and relative to the original code **18.98 s → 4.00 s**, that is **4.75× in total**.
The stage-6 timings were taken in a separate session (see the note at the start of section
5): compare ratios across sections, not absolute seconds.

Memory grew by **1 %**. Through stage 4 the results on deterministic targets were
**bit-identical** to the original code; after stage 6 that no longer holds — a different
summation order in `float` gives differences of order 10⁻⁴ on the weakest significant
matrix elements (section 4, stage 6).

In addition, a **read of uninitialised memory** in the iteration was fixed.

---

## 2. Starting point

The project implements the Discrete Dipole Approximation: light scattering by an
ensemble of dipoles. The core is an iterative solution of the system
`(I − αA)P = αE`, where the product by `A` is evaluated as a three-dimensional
convolution via FFT. The problem size is set by the lattice `i_max × j_max × l_max`;
cubes of 64³ and 128³ are used in practice. Internally the convolution works on an
array zero-padded by a factor of 2 along each axis (that is, 8·N cells).

Call structure:

```
main.cpp          loop over samples / orientations
  particle.cpp    particle shape generation
  dda_func.cpp    DDA_FFT: the iteration
     Convolution  -> FFT_3D / iFFT_3D  (fft_func.cpp)
     LinearComb, Copying, ldcScalMul, ldError   (per-cell loops)
     CabsCext, FarField
```

---

## 3. Where the time actually goes

The profile was taken by instrumenting the current build (64³, 20 steps, 1 sample,
`W=1`, 9.08 s):

| Component | Time | Share |
|---|---|---|
| **FFT** (`FFT_3D` + `iFFT_3D`) | 7.27 s | **80 %** |
| Pointwise multiply in the convolution | 0.61 s | 7 % |
| `CopyingBack` | 0.31 s | 3 % |
| Everything else (per-cell loops, `FarField`, step 0) | ~0.9 s | 10 % |

Inside the FFT: the **1-D kernel is ≈ 63 %**, the line gather/scatter and the
twiddle multiply ≈ 37 %.

Memory: **673 bytes per lattice cell** in total, of which

| Array | Bytes/cell | Comment |
|---|---|---|
| `cmA` (DFT of the interaction kernel) | 384 | shared by all threads |
| `cvWork` (padded convolution buffer) | 192 | per thread |
| `cvP`, `cvP_min`, `cvR`, `cvS` | 96 | per thread |
| `cvG`, `cvH` (solver 2 only) | 48 | per thread |
| `cAlpha` | 1 | |

Hence: 128³ ≈ 1.4 GB with a single thread; 256³ (11.3 GB) does not fit in 8 GB at all.

---

## 4. What was done

### Stage 0 — build (`df8a07a`)

* `compile` switched to `-std=c++17 -pthread`. Apple clang was previously compiling
  the files as **C++98**; without C++17 neither `std::thread` nor lambdas are available.
* 12 cases of `double → float` narrowing in initialiser lists were removed (in
  `dda_func.cpp` and `particle.cpp`) — under C++11+ this is a compile error, and
  **in GCC as well**, so the fix is needed for Intel/AMD too.
* 4 `printf("%d", time_t/clock_t)` warnings were changed to `%ld` (portability).
* `-march=native` is deliberately **not** used: it was measured to give 0 %, and
  leaving it out keeps the binary portable.

### Stage 1 — bit reversal (`a3cf165`, `fft_func.cpp`) — **2.85×**

```c
/* before: malloc + free + pow() for EVERY element */
bits = malloc(number_of_bits*sizeof(short int));
...
(*number_new) += (int)pow(2., (double)i) * bits[i];

/* after */
for (i = 0; i < number_of_bits; i++)
    r = (r << 1) | ((number_old >> i) & 1);
```

The function is called once per FFT point, that is, for a 128³ transform —
**6.3 million times** per three-dimensional transform, and twice per iteration
(forward plus inverse). Microbenchmark: 69 ns/call → 3.5 ns/call. On the real
case: 29.6 s → 10.4 s, output **bit-identical**.

### Stage 2 — memory and a pool of workspaces (`bc67c5c`) — another 1.17×

* `alloc_field3()`: an `ni×nj×nl` array is allocated as **one block** (3 `malloc`
  calls instead of ~33 000 per call), with the row pointers carved out of that
  block. The `cvector***` API is preserved, so no existing loop had to change.
* `dda_workspace` — the seven working arrays live for the whole run instead of
  being allocated and freed twice per sample. This removed repeated page-fault
  storms and was a prerequisite for stage 4.

**An important finding.** Making the storage contiguous on its own produced a
**6 % regression**: the lattice dimensions are powers of two, so every row landed
at a power-of-two offset and the three FFT passes kept hitting the same cache
sets. Adding `FIELD_PAD = 1` (one element of row padding, 24 bytes) turned
−6 % into **+17 %**, at a cost of +0.8 % memory. PAD = 1, 2, 4 and 8 were tried;
one element is enough.

### Fixing a defect (`cc64d5c`)

At the end of `DDA_FFT` there is a fallback to the "best" iterate:

```c
if (Error > Error_min) cvP = cvP_min;
```

but `cvP_min` was filled **only if some iteration improved the discrepancy**. If
the process never improved (typical for the simple iteration method with a large
weight factor, and possible for the other methods on hard targets), the fallback
copied **uninitialised memory** into the dipole moments. The symptom: the same
binary with the same arguments returned `Cabs` = 4.5e7, 2.0e6 and 5.2e10 in three
consecutive runs.

The fix is to store the step-0 iterate immediately after it is computed.
Behaviour is unchanged whenever at least one iteration improves the discrepancy
(in that case `cvP_min` is overwritten anyway).

Consequence for comparing with old results: **the output of solver 3 no longer
matches what the original code produced, and it should not** — it used to be heap
garbage (different from run to run), now it is the step-0 iterate. To verify: for
a diverging iteration the result of `Max_Step=5` is now bit-identical to the
result of `Max_Step=0`, and repeated runs agree with each other. Solvers 0, 1 and
2 are unaffected by the fix (verified bit-identically on 32³, 64³, 128³ and on a
non-cubic 32×64×32 lattice).

### Stage 4 — solving the samples in parallel (`5f81e60`, `ea7530d`, `6b6f799`) — 2.2×

The samples are independent, so the loop over them was parallelised:

* **4a**: `rand()` keeps its state in a global variable and is not thread-safe. A
  thread-local generator (xorshift64*) was added behind `dda_rand()`, returning the
  same 31-bit range; all call sites in `particle.cpp`, `dda_func.cpp` and `main.cpp`
  were switched to it.
* **4b**: an `iVerbose` flag was added to `dda_workspace` — with several threads the
  per-step discrepancy lines are not printed (they would interleave), and one
  summary line per finished sample is printed instead.
* **4c**: the body of the sample loop was moved into a lambda; inside it the local
  names are bound to the worker's private arrays, so **the loop body itself is
  unchanged** — only the accumulators and the printing differ. Workers take sample
  indices from a `std::atomic`; the accumulation of the sums and the update of the
  `.dst` file happen under a mutex once per sample (the file keeps its meaning, so
  a run can still be interrupted and resumed).

Private to a worker: the occupancy map, the polarisabilities, the scattering
geometry, the `Einc`/`Esc` fields, and its `dda_workspace`. Shared read-only:
`cmA`, the FFT tables, `ksc_ini`/`vNormal_ini`.

### Stage 6 — own radix-4/8 1-D FFT kernel (branch `perf/fft-radix48`)

This is the item that section 11 of the previous report listed as the last remaining
lever. The old kernel (`fft_func.cpp`) is a radix-2 DIT that (a) recomputed the bit
reversal on every call and (b) after each stage, starting from the second, copied the
whole line from one array to the other. For a line of N points that is `log2(N)−1` full
copies per transform plus a separate permutation pass — memory traffic, not arithmetic,
was the limit.

What was done (`fft_kernel.hpp`, header-only; the old `FFT`/`iFFT` stay in `fft_func.cpp`
as the fallback path):

* A plan per size N is built once and then only read: a mixed-radix digit-reversal
  permutation for the chosen radix sequence, plus a full table of N `W_N^k` values
  computed directly with `cos`/`sin` in `double` and cast to `float`. Plans are cached
  (up to 16 sizes) and built under a mutex, because `FFT_3D` is called from several
  workers at once.
* The stages run **in place**, in the output buffer: no ping-pong between stages and no
  intermediate `cvWork2`. The input is not modified, but `in` and `out` must not alias:
  the line is loaded already permuted (`out[k]=in[perm[k]]`). The driver still gathers a
  line into `cvWork` and scatters the result from `cvResult`.
* Radix-8 and radix-4 butterflies (radix-8 is built from two radix-4); the three
  components `x, y, z` are handled together. At `i=0` there is no twiddle multiplication;
  inside radix-8 the `±i` twiddles are free and `±(1±i)/√2` costs two multiplications by
  the constant 1/√2 per component; all other elements take ordinary complex products
  with the tabulated twiddles.
* The `log2(N)` split is chosen from a table (`8·8`, `8·4`, `4·4`, `8·4·4`, …), with the
  large radix in the first, most local stages; radix-2 remains only for N=2.
* The kernel is enabled only for powers of two: for `N = 2^k` the driver takes a plan,
  otherwise (or in a `-DFFT_KERNEL_LEGACY` build) the old radix-2 kernel runs. The inverse
  transform uses the same kernel as `conj → DFT → conj`, with no `1/N` factor, exactly like
  the old `iFFT`.

The contract is unchanged: the kernel still computes **`DFT_N` in natural order**. The
"half-spectrum" convention lives in `FFT_3D`, not in the kernel: the padded line of
length `2N` has only `N` non-zero samples, so `X[2l] = DFT_N(x)[l]` and
`X[2l+1] = DFT_N(x·W_{2N}^j)[l]`, and the driver interleaves even/odd bins. A delta test
confirms this: the old kernel is a plain `DFT_N`, while the composite line transform is a
full `DFT_2N` of the zero-padded line.

Measurements (same machine as the rest of the report):

| What | Before | After | Speed-up |
|---|---|---|---|
| 1-D kernel, N=16 | 131 ns/line | 77 | 1.71× |
| 1-D kernel, N=32 | 282 | 184 | 1.53× |
| 1-D kernel, N=64 | 632 | 390 | 1.62× |
| 1-D kernel, N=128 | 1459 | 940 | 1.55× |
| 1-D kernel, N=256 | 3299 | 2121 | 1.56× |
| 64³, 20 steps, 1 sample, `W=1` | 5.67 s | **4.00 s** | **1.42×** |
| 128³, 5 steps, 1 sample, `W=1` | 16.75 s | **11.56 s** | **1.45×** |
| 64³, 4 samples × 5 steps, `W=4` | 2.23 s | **1.58 s** | **1.41×** |
| 64³, original code → this branch | 18.98 s | **4.00 s** | **4.75×** |

The absolute nanoseconds and the ratios in the table are from the stage-6 session. A repeat
run of `compare` on the same machine gives 1.5–1.6× for N=16…256; the absolute values vary
by a few per cent with the machine load.

The estimate from the old section 11 (10–15 % of the end-to-end time) was too
conservative: the kernel got 1.5–1.7× rather than 1.2–1.5×, which at an ≈80 % FFT share
gives 1.42×. That is exactly the Amdahl limit for a 1.6× kernel:
`0.8/1.6 + 0.2 = 0.70`.

**Accuracy.** The kernel itself is more accurate than the old one (error against a
`double` DFT at N=128: 1.1·10⁻⁶ → 8.9·10⁻⁸), because the twiddle table is computed
directly instead of by recurrent multiplication. But the result is no longer
**bit-identical**: `real` is `float` here, and a different summation order gives different
rounding noise. Verified:

* delta test of the new kernel at every position, N=4…256: error ≤ 7.3·10⁻⁸;
* `float` master vs `double` master and `float` new vs `double` master at 64³: the worst
  significant elements differ by 1.2·10⁻⁴ and 1.9·10⁻⁴ — both `float` binaries are equally
  far from the exact answer, and 8.8·10⁻⁵ from each other; the difference grows with the
  number of iterations (1.2·10⁻⁶ for 1 step → 5·10⁻⁵ for 10), as rounding accumulation
  should;
* a build with `-DFFT_KERNEL_LEGACY` (old kernel) reproduces master **bit for bit** — the
  regression test that the driver rework changed nothing.

The checks and measurements are reproduced by the `bench/` harness (`sh bench/build.sh`;
modes `delta`, `newdelta`, `accuracy`, `compare`, `half`, `speed`, `roundtrip`), and the
`.dat` output of different builds is compared with `bench/cmp_dat.py`. It is not part of
the main binary: `compile` globs `*.cpp` in the root, while the harness lives in a
subdirectory.

---

## 5. Results

*Note.* The tables in this section are the stage 0–4 measurements, taken in one session.
The stage-6 timings (section 4) were taken later, in a separate session and on a less busy
machine, so absolute seconds are not comparable across sections — only ratios are. In this
session `master` takes 5.67 s on 64³/20 steps (not 8.85 s as in the table below), and the
new stage 6 takes 4.00 s, i.e. the same 1.42×.

### End-to-end measurements

| Configuration | Before | After | Speed-up |
|---|---|---|---|
| 32³, 19 angles, 40 steps | 6.33 s | 1.83 s | 3.46× |
| 64³, 19 angles, 20 steps | 29.6 s | 8.85 s | 3.34× |
| 64³, 181 angles, 3 steps | 12.2 s | 4.03 s | 3.02× |
| 128³, 19 angles, 3 steps, solver 2 | 120.3 s | 35.7 s | 3.37× |
| 64³, 20 samples × 20 steps, `W=1` | 576.4 s | 171.1 s | 3.37× |
| 64³, 20 samples × 20 steps, **`W=7`** | 576.4 s | **76.6 s** | **7.5×** |

### Scaling with the number of workers (64³, 20 samples × 20 steps)

| W | 1 | 2 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| time, s | 171.1 | 120.9 | 100.3 | 88.7 | 84.6 | **76.6** | 80.4 |
| speed-up | 1.00 | 1.42 | 1.71 | 1.93 | 2.02 | **2.23** | 2.13 |

128³ (4 samples, 2 steps): `W=1` 105.2 s, `W=2` 74.5 s (**1.41×**), `W=4` 78.7 s.

**Why the optimum on 128³ is 2 workers.** A worker's working set there is 674 MB;
two workers already saturate the M1 memory bandwidth (~68 GB/s), and the third and
fourth do not help. On 64³ the set is three times smaller, so 7 workers pay off
(the 8th competes with the system).

**Threads cost nothing:** 4 threads complete the same work in exactly the time
4 independent single-threaded processes take (100.3 s versus 100.5 s).

### Memory

Peak RSS on 128³: 1.508 → **1.524 GB** (+1 %). The increase is only the row padding.

---

## 6. Measured and rejected (with the evidence)

This section matters as much as what was done: it records where not to go back to.

### 6.1. Compacting the fields to active cells — capacity only, not speed

A sphere fills 52 % of its bounding box, so storing `cvP/cvR/cvS/cvP_min` only for
active cells would save −16 % of memory per worker. But on 64³/128³ memory **is not
the constraint** (1.4 GB out of 8 GB), and it would not unlock a 256³ lattice
anyway: `cmA` alone is 6.4 GB there. There is no speed gain either: the per-cell
loops live in the ~10 % of the time that is left after the FFT.

### 6.2. `cmA` cannot be reduced — verified numerically

The hypothesis: the kernel has the form `A(r) = f(|r|)·I + g(|r|)·dyad(r)`, so its
Fourier transform ought to be `Â(k) = A(k)·I + B(k)·k̂k̂` — only 2 complex numbers
instead of 6 (48 → 16 bytes per cell, −36 % memory and −67 % of the convolution
traffic). A diagnostic program was written against the real `CreateAMatrix` +
`FFT_3D`:

```
grid 16³: max relative deviation from A*I+B*kk^ = 2.19e+01
          (only k inside |k|d0 < pi/2)         = 8.67e-01
grid 32³: max relative deviation                = 1.38e+01
          (only k inside |k|d0 < pi/2)         = 1.04e+00
```

The deviation is ~100 % even at the centre of the Brillouin zone: the cubic
anisotropy of the lattice-sampled kernel is comparable to the isotropic part. The
identity **does not hold**, the idea is closed.

### 6.3. Not zeroing the unused part of `cvWork` — 0 %

For FFT mode 0 the transform really does not read 7/8 of the padded buffer
(checked against the loop bounds of all three passes). Removing the zeroing gave
**8.85 s versus 8.81 s** with a bit-identical result. The zeroing loop is cheap and
cache-friendly as it is.

### 6.4. Stage 5, items 1–2: SoA and batched lines — 0.37×…1.01×

Six variants of the 1-D kernel were built and measured separately:

| Variant | ns/line | Speed-up |
|---|---|---|
| A — current code (AoS `cvector`) | 1209 | 1.00× |
| B — SoA by plane (6 arrays) | 2581 | **0.37×** |
| C — AoS, the 6 consecutive `float`s of a butterfly | 1220 | 1.00× |
| E — integer division `i*n/iStep` replaced by a multiply | 1193 | 1.01× |
| all of the above with `-march=native` | — | 1.00× |
| D — batch of 4 lines, layout `[l][4]`, vector along lines | 1132 | **0.90×** |
| D + forced vectorisation (`#pragma clang loop vectorize`) | 1130 | 0.91× |

The reasons: with a 64-element line, spreading the data over planes destroys
locality (a butterfly touches 6 locations instead of two adjacent blocks); and in
most stages the inner loop has 1–2 iterations while the parallelism lies along the
line index, where the data is strided — the compiler states outright
`cost-model indicates that vectorization is not beneficial`.

### 6.5. pocketfft — not applicable (because of the problem structure, not quality)

The library was taken from `github.com/mreineck/pocketfft` (BSD-3, header-only,
portable).

*A note on the integrity of the file*: lines such as `cb.i=twai*t4.r twbi*t3.r;`
look like corruption, but this is the author's idiom — the sign comes from the
macro argument (`POCKETFFT_PARTSTEP5a(...,+tw1i,+tw2i)`). The file is intact, it
compiles, and it computes correctly for all sizes (verified for n = 2, 4, 5, 7, 8,
10, 11).

Measurements:

| Integration | 128³, 3 components |
|---|---|
| our `FFT_3D` (as is) | **71 ms** |
| pocketfft `c2c` over the whole padded volume | 102 ms |

pocketfft is **1.4 times slower**, because our code knows about the zero padding:
on the first pass it reads only the meaningful block (64³ out of 128³) and expands
the result, that is, it does roughly **four times less work**. A general-purpose
library has to transform the whole volume.

Replacing only the inner kernel is not possible either: a delta test shows that
the "1-D FFT" in this code is an N-point transform, and the full spectrum of the
zero-padded line of length 2N is assembled from two of them (hence the `2*l` /
`2*l+1` interleaving in the caller):

```
delta at position 1:
  current kernel: +0.92-0.38i  +0.71-0.71i  +0.38-0.92i     step pi/8  -> 2n-point DFT
  pocketfft:      +0.71-0.71i  +0.00-1.00i  -0.71-0.71i     step pi/4  -> n-point DFT
```

A correct replacement would have to compute the 2n-point transform — twice the
work, which cancels the kernel's advantage (pocketfft is faster "per unit of work"
by only 1.17× at n=64 and 1.48× at n=128).

**Caution:** during this check a prototype with the swapped kernel showed "+3 % on
64³ and +19 % on 128³", but those numbers are **invalid** — it was computing a
different transform. The delta test is what exposed it. Do not use those numbers.

The library header (`pocketfft_hdronly.h`, BSD-3, 121769 bytes, sha256
`518a95c5…44c86d`) is not part of the repository and is not distributed with the
project — it was needed only for this check. To repeat the measurements, take it
from `github.com/mreineck/pocketfft`, branch `cpp`.

---

## 7. Defects found

### 7.1. Uninitialised read in `cvP_min` — **fixed** (`cc64d5c`)

See section 4. Symptom: a non-deterministic result from the simple iteration
method (`Cabs` = 4.5e7 / 2.0e6 / 5.2e10 in three runs).

### 7.2. `iMode_FFT = 1` (Rader–Brenner) is asymptotically poor — **not fixed**

Inside the 1-D `FFT` this mode has a double loop `for(i < N/2) for(j < N/2)`, that
is, **O(L²) per line** and O(N²) per three-dimensional transform. It is already
slower than mode 0 on a 32³ lattice (1.36 s versus 1.00 s), and the gap grows with
size.

**Recommendation: do not use `iMode_FFT = 1`.** Mode 0 (Cooley–Tukey) is the only
practical choice. The code of mode 1 is left untouched so that behaviour does not
change.

### 7.3. Building as C++98 by default — **fixed** (`df8a07a`)

### 7.4. `DDA_FFT` was not reentrant — consequence fixed along the way

Before stage 2 the working arrays were allocated and freed inside `DDA_FFT`, which
by itself ruled out parallelism. The state now lives in `dda_workspace`.

### 7.5. `FourierPhase(1)` writes past the end of the array — **not fixed**

Found while debugging the stage-6 harness. `FourierPhase(i_max, cW)` unconditionally
fills `cW[0]` and `cW[1]`, while the `cW_2_*` tables are allocated with `i_max/2`
elements. For a lattice dimension of 2 that is `FourierPhase(1)` into a one-element
array — an 8-byte write past the end of the heap block. It is reachable only for a
lattice with dimension 2 (physically meaningless), and only in the old code; the new
kernel does not use that table, but `main.cpp` still builds it. Left untouched so as
not to change behaviour beyond the task.

### 7.6. Particle generation can hang — **not fixed**

Found while regressing the new interface: all 14 target types were compared with
the reference build byte for byte, and on two parameter sets **both builds hang** —
the old one (positional arguments) and the new one. So the defect pre-dates the
command-line work and is unrelated to it.

**`iGenerateHedgehog` (`particle.cpp`)** — the same placement scheme as in
`iGenerateAnthrax` and `iGenerateDeadAnthrax`:

```c
while(m<iWormsNumber)
{   cl0.r=-1.;
    while(cl0.r<0.)      /* the cylinder does not fit — draw another one */
    {   ...
        cl0=worm(...,length_worm,radius_worm,R2);
    }
    for(j=0;j<m;j++) { l=CrossTwoWorms(cl0,worms[j]); if(l==1) break; }
    if(l==0) { worms[m]=cl0; m++; }
}
```

Both loops are unbounded when no valid placement exists: the inner one when the
monomer (`--monomer-radius` together with `--monomer-length-ratio`) does not fit
into the matrix, the outer one when `iWormsNumber`, derived from `--packing`, is
larger than the number of non-overlapping worms that fit. There is no attempt
limit, so the call never returns.

**`iGenerateHierarchicalEllipsoid`** — the seed-cell search:

```c
do
{   iFlag=0;
    do
    {   v0 = random cell of the lattice;
        i0 = cAlpha[(int)v0.x][(int)v0.y][(int)v0.z];
    }
    while(i0!=1);        /* a cell of the first material is required */
    ...  /* reject the already chosen seeds */
}
while(iFlag==1);
vSeeds[m]=v0;
```

The inner loop never ends when the current level has **no cell with value 1**
(or fewer than `iNoSeeds` such cells). That is typical for a tight matrix: for an
ellipsoid with small semi-axes on an 8³ lattice the inner region is empty.

**Reproduction** (the old build hangs as well with the same parameters):

```sh
# type 7: packing too high for an 8³ lattice
./dda_fft.out --target hedgehog --packing 0.3 --monomer-radius 1 \
    --monomer-length-ratio 1 --material 1.5 0.01 --grid 8 --size-parameter 2

# type 13: the 0.8 x 0.6 ellipsoid leaves no first-material cells on 8³
./dda_fft.out --target hierarchical --aspect 0.8 0.6 --levels 2 \
    --material 1.5 0.01 --material 1.33 0.0 --seeds 111 222 \
    --grid 8 --size-parameter 2
```

With parameters that leave enough room both types run normally: type 7 with
`--packing 0.05 --monomer-radius 0.5 --grid 16`, type 13 with
`--aspect 1 1 --levels 2 --grid 16` (the new and old builds agreed byte for byte
on those sets).

**Why it is not fixed.** The proper fix is to bound the number of attempts and
return an error or a partially filled particle. That changes the geometry and
therefore the numbers, so it needs a separate task with a physics check. Until
then: choose `--packing`/`--monomer-*` that fit the matrix volume, use a large
enough lattice for type 13, and put an external timeout around scripted runs.

---

## 8. How to use

### Build

```sh
sh compile
```

A C++17 compiler is required: GCC 7+ or clang 5+. On Linux/Intel/AMD the system
`g++` is enough; no extra packages are needed: argtable3 is vendored in
`third_party/argtable3/` (BSD-3-Clause) and built as C by the same `compile` script.

### Number of threads

```sh
./dda_fft.out --threads 4 <options>            # explicit
./dda_fft.out <options>                        # default: one per core
```

The `DDA_THREADS` environment variable is still honoured when `--threads` is
not given.

By default the core count is used, capped by the number of samples and by what
memory allows (70 % of physical RAM minus the shared `A` matrix). The choice is
reported as

```
 Worker threads: N
```

With `N = 1` the output is exactly as before (per-step discrepancies); with
`N > 1` there is one summary line per finished sample.

**The best number of threads depends on the lattice size**, because the limit is
memory bandwidth, not cores. On an M1 with 8 GB the measured optimum was: 64³ → 7
threads, 128³ → 2 threads. On other machines the optimum differs; one run with 1,
2, 4 and 8 threads is enough to find it.

### Running

```sh
# 64³ sphere, m = 1.5 + 0.01i, 19 angles, 20 steps
./dda_fft.out --threads 7 --target sphere --material 1.5 0.01 --grid 64 \
              --size-parameter 20
```

The FFT method and the solver are selected with `--fft` and `--solver`.

---

## 9. Portability (Intel / AMD / Apple Silicon)

* The changes are architecture-independent: no intrinsics and no `-march=native`.
  The syntax was checked for x86-64
  (`clang++ --target=x86_64-apple-darwin -fsyntax-only`).
* `FIELD_PAD = 1` also cures cache-set conflicts on x86 (4 KB pages, 64-byte
  lines): the row stride stops being a multiple of a power of two.
* **`long double` is not a problem on x86.** On Apple Silicon it equals `double`
  (8 bytes); on x86 it is 80 bits (x87, no vectorisation). The share of the three
  accumulation functions (`ldcScalMul`, `ldScalMulConj`, `ldError`) was measured:
  **0.075 s out of 9.5 s = 0.8 %**. Even a 3–4× slowdown of x87 means +2–3 % on
  Intel. Not worth touching.
* **Results will not be bit-identical between arm64 and x86** — because of the
  different width of `long double` in the accumulators. This was already true of
  the code; within one platform everything is deterministic.

---

## 10. Reproducing the measurements

General form of the command:

```sh
./dda_fft.out --state out.dst --output out.dat --samples <N> \
    --orientation <fixed|mc> --target sphere --material <mre> <mi> \
    --grid <i> <j> <l> --size-parameter <x> --polarizability <1..5> \
    --tolerance <Error_max> --max-steps <Max_Step> --fft <0|1> --solver <0..3>
```

The measurements quoted in this report (19 angles; sphere; `m = 1.5+0.01i`;
`x = 20`; polarisability 1):

```sh
# 64³, 20 steps, solver 0
./dda_fft.out --state o.dst --output o.dat --target sphere --material 1.5 0.01 \
    --grid 64 64 64 --size-parameter 20 --tolerance 1e-3 --max-steps 20 --fft 0 --solver 0

# 128³, 3 steps, solver 2
./dda_fft.out --state o.dst --output o.dat --target sphere --material 1.5 0.01 \
    --grid 128 --size-parameter 40 --tolerance 1e-3 --max-steps 3 --fft 0 --solver 2

# 64³, 20 samples × 20 steps (the main "user" scenario)
./dda_fft.out --threads 7 --state o.dst --output o.dat --samples 20 \
    --target sphere --material 1.5 0.01 --grid 64 --size-parameter 20 \
    --tolerance 1e-3 --max-steps 20 --fft 0 --solver 0
```

The correctness checks used during the work:

```sh
# bit-identical result
cmp new.dat old.dat

# thread safety
cc -O1 -std=c11 -c third_party/argtable3/argtable3.c -o argtable3.o
g++ -O1 -g -std=c++17 -pthread -fsanitize=thread  -o t *.cpp argtable3.o && ./t --threads 4 ...
g++ -O1 -g -std=c++17 -pthread -fsanitize=address -o a *.cpp argtable3.o && ./a --threads 4 ...
```

Checked: solvers 0/1/2 (and 3 — see above) and both FFT modes; lattices 32³, 64³,
128³ and the non-cubic 32×64×32; fixed and random orientations; ThreadSanitizer
and AddressSanitizer — clean.

---

## 11. What is left and why it was not done

The radix-4/8 kernel — the only significant lever that was left — is closed by stage 6.
After it the FFT convolution hits the Amdahl limit: with an ≈80 % FFT share and a 1.6×
kernel, the other 20 % (the pointwise multiplication, the row gather/scatter,
`CopyingBack`, step 0) cap the gain at ≈1.43×, and that is already reached.

Small reserves — fusing the row gather with the permutation, moving the line buffers into
`dda_workspace` instead of `malloc` per 3-D call — are worth a few per cent and do not pay
for the extra complexity.

Further gains are only possible by changing the algorithm (for instance multi-level or
hierarchical methods instead of the FFT convolution on large lattices), not by optimising
the present scheme.

---

## 12. Change history

```
af63403  Baseline: DDA/FFT source as received
df8a07a  Stage 0: build with C++17, fix narrowing inits and printf format specifiers
a3cf165  Stage 1: replace malloc/pow bit-reversal with shift-based version
bc67c5c  Stage 2: contiguous field storage and a reusable solver workspace
cc64d5c  Fix: keep the step-0 iterate in cvP_min
5f81e60  Stage 4a: per-thread random number generator
ea7530d  Stage 4b: optional quiet mode for the solver
6b6f799  Stage 4c: solve the samples in parallel
(perf/fft-radix48)  Stage 6: mixed-radix 4/8 in-place 1-D FFT kernel
4eae0b7  Vendor argtable3 3.3.1 for command-line parsing
a0f77ac  Replace positional arguments with named command-line options
```

In total relative to the base: **8 files, +688 / −388 lines** for stages 0–4; stage 6 adds
`fft_kernel.hpp` and a change to `fft_func.cpp` (+67 / −46). The changes do not touch the
physics — only the build, the memory layout, the way the FFT is computed, and the order in
which independent samples are traversed.

The program banner in `main.cpp` is updated to version **7.1**, naming the person responsible
for the modification and optimisation (`K.A.Shmirko`, <kshmirko@dvo.ru>, 2026).

The `feature/cli-argtable3` branch (tag **v7.2**) moves the configuration to named
command-line options: the parsing lives in `options.cpp`/`options.hpp`, the positional
arguments are gone, and range checks, `--seed` and `--restart` were added. The version in
the banner, in `DDA_VERSION` and in the documentation is **7.2**.

The only third-party code in the repository is the vendored argtable3 v3.3.1
(`third_party/argtable3/`, BSD-3-Clause); the pocketfft header mentioned in section 6.5 is
not part of the distribution.
