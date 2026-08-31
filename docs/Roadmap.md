# Plan de paralelizacion y benchmarking

Este documento es el plan de trabajo para ir de la version secuencial a la
paralelizada, midiendo speedup en cada paso. La idea es ir llenando las
tablas de resultados a medida que se prueba cada version/rama.

## Mapa de ramas

| Rama            | Binario           | Hilos | Solver de presion/difusion      | Fullscreen | CLI / codigo |
|-----------------|-------------------|-------|----------------------------------|------------|--------------|
| `master`        | `Screensaver_seq` | 1     | Gauss-Seidel fila por fila       | si (`-p`)  | vieja, sin tocar (congelada en `f6bca8a`) |
| `sequential`    | `ss`              | 1     | Gauss-Seidel fila por fila       | si (`-p`)  | actual: ingles, sin `-t`/`-b`, `src/` plano |
| `parallel-omp`  | `ss`              | N     | Red-black (paralelizable)        | no         | vieja, congelada antes de la traduccion/limpieza |

`master` quedo intencionalmente sin tocar despues del reset a `f6bca8a` (ver
el resto de esta conversacion): todo el trabajo de traduccion a ingles,
aplanado de `src/secuencial/` a `src/`, renombrado del binario a `ss` y
limpieza de la CLI (quitar `-t`/`-b`) paso unicamente en `sequential`. Si se
compara con `master` hay que tener en cuenta que ese binario todavia se
llama `Screensaver_seq`, sigue en espanol y sigue teniendo `-t`/`-b`.

`sequential` = referencia correcta para medir speedup: mismos parametros,
mismo comportamiento visual, sin OpenMP y con el algoritmo de solver que de
verdad es secuencial (no el red-black, que solo existe para poder
paralelizarse — ver la nota en `solver.c`).

`parallel-omp` es la version con OpenMP tal cual quedo tras el refactor de
render/fuentes/fondo, con `-fopenmp` en el Makefile. **A partir de aqui se
deja como referencia congelada** (para ver el resultado final de a donde se
quiere llegar) y todo el trabajo de optimizacion incremental se hace sobre
`sequential`, creando una rama nueva por cada paso del checklist mas abajo.
`parallel-omp` tambien quedo congelada con la CLI vieja (todavia tiene `-t`
para dt y `-b` para alternar n-body, ambos ya eliminados en `sequential`: ahi
dt es fijo y el n-body siempre esta activo) y sin fullscreen — no se le sigue
actualizando a la par de `sequential`. Al comparar FPS entre ramas, tener en
cuenta que la CLI difiere.

**Ninguna rama de paso incremental existe todavia.** Todo el trabajo activo
pasa en `sequential`: se implementa el paso ahi, se corre el benchmark, se
llena la tabla, se hace commit, y **recien entonces** se crea la rama de
ese paso como snapshot (ver "Roadmap de pasos" mas abajo para el orden y
nombres planeados: `01-rb-tree`, `02-omp-loops`, `03-omp-solver`,
`04-schedule-tuning`, `05-collapse-tuning`). No crear ramas por adelantado.

---

## Parametros (linea de comandos y constantes)

### Flags de `./ss`

| Flag | Rango / default | Efecto |
|------|------------------|--------|
| `-n <N>` | `[16..1024]`, def. `256` | Resolucion de la malla (N x N celdas interiores). El costo del solver escala ~O(N^2) por iteracion; es el parametro mas importante para medir speedup. |
| `-f <F>` | `[1..256]`, def. `6` | Cantidad de fuentes de tinta. |
| `-W <ancho>` | min. `640`, def. `1920` | Ancho de ventana en px (no afecta el costo del solver, solo el de render). |
| `-H <alto>` | min. `480`, def. `1080` | Alto de ventana en px. |
| `-s <semilla>` | def. reloj del sistema | Semilla PRNG — fijarla para comparaciones reproducibles entre corridas. |
| `-v <visc>` | `[0.0..1.0]`, def. `0.0` | Viscosidad cinematica. |
| `-d <diff>` | `[0.0..1.0]`, def. `0.0001` | Coeficiente de difusion de tinta (el "blur"/suavizado). |
| `-p`, `-F` | desactivado por default | Pantalla completa (solo en `sequential`/`master` por ahora). |
| `-h` | — | Ayuda. |

`dt` ya no es ajustable por CLI (queda fijo en `DT_DEFAULT`) y el sistema de
n-cuerpos siempre esta activo (no hay flag para desactivarlo) — ambos
simplificados fuera de la rama `parallel-omp`, que sigue teniendo `-t` y
`-b` por ser una rama congelada con la CLI vieja.

### Constantes (`include/common.h`, `include/sources.h`, `src/nbody.c`)

Estas no son ajustables por CLI, pero son relevantes si se quiere variar el
"blur"/comportamiento en experimentos mas finos (requiere recompilar):

| Constante | Valor | Que hace |
|---|---|---|
| `DT_DEFAULT` | `0.07` | Paso de tiempo de la simulacion (ya no ajustable por CLI en `sequential`). |
| `GAUSS_SEIDEL_ITERS` | `20` | Iteraciones del solver lineal por llamada — sube linealmente el costo de `solve_linear()`, que es el hotspot principal. |
| `DISSIPATION` | `0.995` | Retencion de tinta por frame. |
| `BRIGHTNESS_FACTOR` / `CONTRAST_FACTOR` / `SATURATION_FACTOR` / `ALPHA_GAMMA` | `0.6` / `1.0` / `2.4` / `0.7` | Solo afectan el render (color), no el costo computacional. |
| `SOURCE_RADIUS_BASE` / `SOURCE_SIGMA_BASE` | `2` / `1.1` | Tamano del vecindario de inyeccion gaussiana por fuente; escala con `-n` via `source_scale()`. |
| `NBODY_G` / `NBODY_SOFTENING` / `NBODY_VEL_MAX` / `NBODY_RESTITUTION` | `0.3` / `6.0` / `5.0` / `0.9` | Fisica del sistema de n-cuerpos (siempre activo); el costo es O(F^2) en cantidad de fuentes. |

**Parametros clave para el plan de benchmarking**: `-n` (resolucion de malla,
el que mas importa), `-f` (cantidad de fuentes, que ahora siempre mueve el
costo O(F^2) del n-body ya que esta siempre activo), y el numero de hilos
OpenMP (`OMP_NUM_THREADS`, ver abajo — no es un flag del programa, es una
variable de entorno).

---

## Catalogo de optimizaciones (parallel-omp)

Cada fila es un sitio de paralelizacion real en el codigo. "Como
implementarlo" describe el patron para poder reproducirlo si se hace una
rama incremental que las va agregando una por una.

### 1. `add_source()` — `src/sources.c`

- **Que hace**: suma el termino fuente (velocidad/tinta inyectada) a cada
  celda, `O(N^2)`, sin dependencias entre celdas.
- **Como se paralelizo**: `#pragma omp parallel for schedule(dynamic, 16)`
  sobre el loop plano de `celdas_total`.
- **Por que funciona sin cambios de algoritmo**: cada celda se escribe una
  sola vez y no depende de otras celdas — paralelizable "tal cual".
- **Como implementarlo**: agregar el pragma justo antes del `for`. Nada mas
  cambia.

### 2. `dissipate_ink()` — `src/sources.c`

- Mismo patron que (1): `#pragma omp parallel for schedule(dynamic, 16)`
  sobre un loop plano e independiente por celda (`tinta *= DISSIPATION`).

### 3. Solver de presion/difusion: Gauss-Seidel -> red-black — `src/solver.c`

- **Que hace**: resuelve el sistema lineal disperso de la difusion implicita
  y de la proyeccion de Hodge, iterativamente (`GAUSS_SEIDEL_ITERS = 20`
  iteraciones por llamada, se llama varias veces por frame). Es el
  **hotspot principal** del programa.
- **El problema**: Gauss-Seidel fila por fila lee valores ya actualizados en
  la misma iteracion (`campo[i-1,j]` ya fue escrito antes de leer
  `campo[i,j]`) — eso crea una dependencia de datos secuencial estricta,
  no paralelizable directamente.
- **La solucion (algoritmica, no solo `#pragma`)**: recorrido *red-black*.
  Se divide la malla en un tablero de ajedrez por paridad `(i+j) % 2`. Los 4
  vecinos de una celda de un color son siempre del color opuesto, asi que
  todas las celdas del mismo color se pueden actualizar en paralelo sin
  pisarse. Cada iteracion pasa a ser: actualizar todas las "rojas", barrera,
  actualizar todas las "negras", barrera.
- **Como implementarlo**:
  1. Separar el cuerpo del loop en una funcion `relax_color(campo,
     campo_previo, a, inverso_c, paridad)` que recorre solo las celdas de
     una paridad (`inicio = ((1+j)%2 == paridad) ? 1 : 2`, luego `i += 2`).
  2. Envolver las `GAUSS_SEIDEL_ITERS` iteraciones en **una sola** region
     `#pragma omp parallel` (abrir/cerrar el team de hilos tiene costo fijo;
     hacerlo una vez en vez de 2*20 veces importa en mallas chicas).
  3. Dentro, cada `relax_color()` usa `#pragma omp for schedule(dynamic,
     16)` (barrera implicita automatica al final).
  4. Aplicar la condicion de frontera con `#pragma omp single` (un solo hilo
     la aplica mientras los demas esperan en la barrera implicita del
     `single`).
- **Nota de correctitud numerica**: red-black converge a un resultado
  practicamente identico a Gauss-Seidel secuencial (de hecho a veces
  converge un poco mas rapido), asi que el comportamiento visual no deberia
  notarse distinto entre `sequential` y `parallel-omp`.

### 4. `advect()` — `src/solver.c`

- **Que hace**: adveccion semi-Lagrangiana, cada celda destino solo *lee*
  `campo_previo` (no lo modifica) e interpola — sin dependencias entre
  celdas de salida.
- **Como se paralelizo**: `#pragma omp parallel for collapse(2)
  private(...) schedule(dynamic, 16)` sobre los dos loops anidados (i, j).
  `collapse(2)` funde ambas dimensiones en un solo rango repartible, para
  que mallas angostas (menos filas que nucleos) sigan repartiendo trabajo
  entre todos los hilos.
- **Como implementarlo**: agregar el pragma con `private()` listando todas
  las variables escalares temporales del cuerpo del loop (`i0,i1,j0,j1,
  origen_x,origen_y,pesos...`) para que cada hilo tenga su propia copia.

### 5. Proyeccion de Hodge (divergencia + gradiente) — `src/solver.c`

- Dos loops mas con el mismo patron que (4): cada celda destino es
  independiente. `#pragma omp parallel for collapse(2) schedule(dynamic,
  16)` en cada uno.

### 6. `render_ink()` — `src/render.c`

- **Que hace**: vuelca la densidad de tinta a la textura SDL, un pixel de
  textura por pixel de pantalla — cada pixel de salida es independiente.
- **Como se paralelizo**: `#pragma omp parallel for collapse(2)
  private(fx, fy, valor_r, valor_g, valor_b, promedio) schedule(dynamic,
  16)`.
- **Detalle de implementacion importante**: para que `collapse(2)` reparta
  trabajo entre hilos incluso en mallas angostas, el calculo de `fila`
  (puntero a la fila de pixeles) y `fy` (que en realidad solo depende de
  `j`, no de `i`) se recalculan dentro del loop interior en vez de estar
  izados fuera del loop de `j` — es un recalculo barato (una resta y una
  multiplicacion) a cambio de mejor paralelismo. **Esta es la unica
  optimizacion de la lista que sacrifica algo de trabajo redundante para
  ganar paralelismo** — en la version `sequential` esto se revirtio a la
  forma hoisted original porque en un solo hilo es puro desperdicio.

### `schedule(dynamic, 16)` — por que se uso en todos los sitios

Se eligio `dynamic` (en vez de `static`, el default) con chunk de 16 filas.
`dynamic` reparte bloques de trabajo bajo demanda segun cada hilo termina su
bloque anterior, en vez de dividir el rango en partes iguales al inicio. Con
`static` un hilo que le toque una zona "mas cara" (mas rara en este
problema, pero puede pasar por variaciones de cache/contencion de memoria
entre iteraciones) deja al resto esperando; con `dynamic` los hilos ociosos
roban trabajo. El chunk de 16 evita el overhead de pedir trabajo celda por
celda.

---

## La estrategia: un paso a la vez, misma bateria de pruebas cada vez

El plan real es: cada optimizacion (algoritmica o de paralelizacion) va en
su propia rama, arrancando siempre desde la ultima rama "buena" (no desde
`sequential` cada vez, sino acumulando). Despues de CADA paso se corre
exactamente la misma bateria de benchmarks (misma metodologia, mismos
valores de `-n`/`-f`/hilos), y se le agrega una columna a la tabla maestra
de resultados. Asi cada numero es comparable con el anterior y se puede
hablar de "de X a Y, +Z% FPS" paso por paso, no solo secuencial-vs-paralelo-
al-final.

El eje que mas importa variar en cada bateria es **`-n`** (tamano de
celda/resolucion de la malla) porque domina el costo (`O(N^2)` por
iteracion, y son `GAUSS_SEIDEL_ITERS` iteraciones por paso, varias veces por
frame). El eje secundario es **`-f`** (cantidad de fuentes/cuerpos): el
n-body siempre esta activo, asi que `-f` mueve directamente su costo
(`O(F^2)` sin arbol, `O(F log F)` con arbol) ademas del costo de inyeccion.
Los demas flags (viscosidad, difusion, ventana, etc.) no cambian el costo
computacional, solo el aspecto visual — no hace falta variarlos en el
benchmark de rendimiento.

### Roadmap de pasos

Orden planeado (nombres de rama sugeridos para cuando se llegue a cada
uno — **ninguna de estas ramas existe todavia**, se crean sobre la marcha):

```
sequential ──▶ 01-rb-tree ──▶ 02-omp-loops ──▶ 03-omp-solver ──▶ 04-schedule-tuning ──▶ 05-collapse-tuning
```

- [x] **Paso 0, baseline** — 1 hilo, n-body O(F^2) con fuerza bruta, solver
      Gauss-Seidel fila por fila. *(ya no es la punta de `sequential`, que
      avanzo al Paso 1 -- este estado sigue vivo en `master`, que se quedo
      congelada ahi a proposito. Columna "baseline" en la tabla maestra.)*
- [x] **Paso 1 — `01-rb-tree`**: arbol para n-body (Barnes-Hut): reemplazo el
      loop O(F^2) de `update_nbody_sources()` (`nbody.c`) por un quadtree
      con aproximacion por centro de masa, O(F log F). Sigue siendo
      secuencial (sin OpenMP todavia). *(implementado y medido -- FPS
      practicamente identico al baseline en todo el rango probado, incluso
      con f alto: a estas resoluciones el solver de fluidos domina tanto
      que el algoritmo de n-body no se nota todavia. Commit vive en
      `sequential`, con `01-rb-tree` apuntando al mismo commit.)*
- [ ] **Paso 2 — `02-omp-loops`**: paralelizar los loops "faciles": los sitios
      1, 2, 4, 5, 6 del catalogo de arriba (fuente/disipacion/adveccion/
      proyeccion/render) — todos son `#pragma omp parallel for` directo,
      sin cambiar el algoritmo. El solver de presion/difusion se queda
      secuencial todavia (sigue siendo Gauss-Seidel fila por fila).
- [ ] **Paso 3 — `03-omp-solver`**: red-black + paralelizar el solver: el
      cambio algoritmico del catalogo (seccion 3 de arriba) — se espera que
      sea el salto mas grande, es el hotspot principal.
- [ ] **Paso 4 — `04-schedule-tuning`**: probar `static`/`dynamic`/`guided`
      y distintos tamanos de chunk sobre el resultado del paso 3, variando
      `-n` (mallas chicas vs. grandes reaccionan distinto al overhead de
      reparto de trabajo).
- [ ] **Paso 5 — `05-collapse-tuning`**: `collapse(2)` si/no, comparar contra
      paralelizar solo el loop externo, especialmente relevante en `-n`
      chico (menos filas que `nproc`).
- [ ] **Paso 6 (opcional)** — SIMD / autovectorizacion, `-O2` vs `-O3` como
      variable de build.

`parallel-omp` (la rama congelada) ya tiene implementados los pasos 2 y 3
juntos (mas los cambios visuales de resolucion/fondo/render que no son
parte de este roadmap de performance) — sirve como referencia de "a donde
se puede llegar", pero la cadena de arriba se construye de nuevo, paso por
paso, para poder medir cada incremento por separado.

**Flujo de trabajo real, un paso a la vez** (todo el trabajo activo pasa en
`sequential` hasta terminar un paso):

1. Sobre `sequential`, implementar el cambio del paso actual.
2. Compilar (`make`) y correr la bateria de benchmark completa (ver abajo).
3. Pegar los resultados en la tabla maestra de este documento.
4. Commit en `sequential` (codigo + tabla actualizada + checkbox marcado).
5. Recien ahi, crear la rama de este paso apuntando a ese commit
   (`git branch <nombre-del-paso>`, p.ej. `git branch 01-rb-tree` despues de
   terminar el Paso 1) como snapshot/checkpoint — no antes.
6. Seguir trabajando en `sequential` para el siguiente paso.

No crear ramas por adelantado para pasos que todavia no se implementaron.

---

## Metodologia de benchmark (igual para cada paso)

1. **Metrica**: el programa imprime FPS por stdout cada ~0.5s
   (`printf("FPS= %.2f\n", ...)`). Correr cada configuracion un tiempo fijo
   (p.ej. 30s) y promediar, ignorando las primeras ~5 lineas (warm-up /
   cache fria).
2. **Reproducibilidad**: fijar siempre la semilla (`-s 42`).
3. **Comando base** (ajustar `-n`/`-f`/flags segun la fila de la tabla):
   ```bash
   timeout 30s ./ss -n <N> -f <F> -s 42 2>/dev/null \
     | grep "FPS=" | tail -n +6 \
     | awk -F'= ' '{s+=$2; c++} END {printf "FPS promedio: %.2f\n", s/c}'
   ```
   (Sin display disponible: envolver con `xvfb-run`.)
4. **Hilos OpenMP** (solo aplica desde el Paso 2 en adelante): controlar sin
   recompilar con `OMP_NUM_THREADS`:
   ```bash
   OMP_NUM_THREADS=4 ./ss -n 256 -f 6 -s 42
   ```
   Nucleos disponibles en esta maquina: `nproc` → **20**. Para los pasos que
   todavia no tienen OpenMP, correr con 1 hilo nada mas (no aplica variar).
5. **Speedup** siempre respecto al paso anterior en la tabla maestra:
   `speedup_paso_N = FPS_paso_N / FPS_paso_(N-1)` (mismos `-n`/`-f`/flags).
   **Speedup acumulado** = `FPS_paso_N / FPS_paso_0` (contra `sequential`).

### Matriz estandar a correr en cada paso

Repetir esta misma matriz completa despues de cada paso del roadmap y pegar
los resultados en la tabla maestra de abajo (agregar mas filas de matriz
si algun paso lo amerita, p.ej. mas hilos una vez que haya OpenMP).

| # | `-n` | `-f` | Hilos | Por que esta en la matriz |
|---|------|------|-------|------------------------------|
| A6  | 64   | 6    | 1     | malla chica, referencia rapida |
| A12 | 64   | 12   | 1     | idem, con el doble de fuentes/cuerpos |
| B6  | 256  | 6    | 1     | default del programa |
| B12 | 256  | 12   | 1     | idem, con el doble de fuentes/cuerpos |
| C6  | 512  | 6    | 1     | malla grande, domina el solver |
| C12 | 512  | 12   | 1     | idem, con el doble de fuentes/cuerpos |
| H6  | 600  | 6    | 1     | rango donde se observo ~20-30 FPS en pruebas anteriores |
| H12 | 600  | 12   | 1     | idem, con el doble de fuentes/cuerpos |
| I6  | 700  | 6    | 1     | idem, extremo superior de ese mismo rango |
| I12 | 700  | 12   | 1     | idem, con el doble de fuentes/cuerpos |
| D6  | 1024 | 6    | 1     | malla muy grande, limite superior |
| D12 | 1024 | 12   | 1     | idem, con el doble de fuentes/cuerpos |
| E   | 256  | 4    | 1     | pocas fuentes/cuerpos |
| F   | 256  | 64   | 1     | muchas fuentes/cuerpos, resalta el costo O(F^2)/O(F log F) del n-body |
| G   | 256  | 256  | 1     | limite superior de fuentes/cuerpos |

#### Comandos (copiar/pegar tal cual, un paso a la vez)

Parado en la rama del paso que se esta midiendo, con el binario ya
compilado (`make`):

La matriz completa (15 filas) de un tiron, imprime una linea por fila,
listo para pegar en la tabla maestra. **Esto SI corre el benchmark
automaticamente** (compila una vez, corre cada configuracion 30s, promedia
el FPS), pero **no escribe en este archivo**: hay que copiar las lineas
impresas a mano a la tabla de abajo.

```bash
for row in "A6 64 6" "A12 64 12" "B6 256 6" "B12 256 12" "C6 512 6" "C12 512 12" \
           "H6 600 6" "H12 600 12" "I6 700 6" "I12 700 12" "D6 1024 6" "D12 1024 12" \
           "E 256 4" "F 256 64" "G 256 256"; do
  set -- $row
  label=$1; n=$2; f=$3
  timeout 30s ./ss -n "$n" -f "$f" -s 42 2>/dev/null | grep "FPS=" | tail -n +6 \
    | awk -v l="$label" -v n="$n" -v f="$f" -F'= ' \
      '{s+=$2; c++} END {printf "%s  n=%-5s f=%-4s -> %.2f FPS\n", l, n, f, s/c}'
done
```

A partir del Paso 2 (`02-omp-loops` en adelante, cuando ya hay OpenMP), repetir
tambien la fila B variando `OMP_NUM_THREADS` para la curva de
escalabilidad:

```bash
for t in 1 2 4 8 16 20; do
  OMP_NUM_THREADS=$t timeout 30s ./ss -n 256 -f 6 -s 42 2>/dev/null \
    | grep "FPS=" | tail -n +6 \
    | awk -v t="$t" -F'= ' '{s+=$2; c++} END {printf "hilos=%-2s -> %.2f FPS\n", t, s/c}'
done
```

---

## Tabla maestra de resultados

FPS promedio por combinacion de fila-de-matriz x paso. Ir agregando una
columna por cada paso conforme se implementa y se corre la bateria.

| Fila matriz | baseline (`master`, pre-Barnes-Hut) | `sequential` / `01-rb-tree` (Barnes-Hut) | `02-omp-loops` | `03-omp-solver` | `04-schedule-tuning` | `05-collapse-tuning` |
|---|---|---|---|---|---|---|
| A6 (n=64, f=6) | 17.21 | 17.19 | | | | |
| A12 (n=64, f=12) | 17.16 | 17.20 | | | | |
| B6 (n=256, f=6) | 10.52 | 10.51 | | | | |
| B12 (n=256, f=12) | 10.95 | 10.95 | | | | |
| C6 (n=512, f=6) | 3.75 | 3.76 | | | | |
| C12 (n=512, f=12) | 4.30 | 4.30 | | | | |
| H6 (n=600, f=6) | 2.69 | 2.69 | | | | |
| H12 (n=600, f=12) | 3.09 | 3.09 | | | | |
| I6 (n=700, f=6) | 1.90 | 1.91 | | | | |
| I12 (n=700, f=12) | 2.10 | 2.10 | | | | |
| D6 (n=1024, f=6) | 0.92 | 0.92 | | | | |
| D12 (n=1024, f=12) | 1.00 | 0.99 | | | | |
| E (n=256, f=4) | 7.52 | 7.52 | | | | |
| F (n=256, f=64) | 11.39 | 11.39 | | | | |
| G (n=256, f=256) | 11.37 | 11.38 | | | | |

Medido en local (display real, no headless), corridas de ~15s por fila,
promedio ignorando los primeros ~2s de warm-up, `-s 42` en todas.
Barnes-Hut da practicamente el mismo FPS que fuerza bruta en todo el rango
probado (incluyendo f=64/256): a estas resoluciones de malla el solver de
fluidos domina tanto el costo que el algoritmo de n-cuerpos, sea O(F^2) o
O(F log F), no se nota en el FPS total. Sigue siendo la base correcta para
cuando se paralelice el solver (Pasos 2+), donde el resto del frame sera
mas barato relativamente y el n-body podria empezar a pesar mas.

### Curva de escalabilidad por hilos (desde el Paso 2 en adelante)

Usar la fila B (n=256, f=6) como referencia, variando `OMP_NUM_THREADS`:

| Hilos | FPS `02-omp-loops` | FPS `03-omp-solver` | FPS `04-schedule-tuning` | FPS `05-collapse-tuning` |
|-------|------------------|-------------------|-------------------------|--------------------------|
| 1     |                  |                   |                         |                          |
| 2     |                  |                   |                         |                          |
| 4     |                  |                   |                         |                          |
| 8     |                  |                   |                         |                          |
| 16    |                  |                   |                         |                          |
| 20    |                  |                   |                         |                          |

---

## Como correr una comparacion rapida entre ramas

Recorre toda la cadena del roadmap (las que ya tengan commits mas alla del
punto donde se crearon daran resultados; las que siguen vacias daran el
mismo numero que su padre, eso es esperado):

```bash
for branch in sequential 01-rb-tree 02-omp-loops 03-omp-solver 04-schedule-tuning 05-collapse-tuning parallel-omp; do
  git checkout "$branch" >/dev/null 2>&1 && make clean >/dev/null && make >/dev/null 2>&1 && \
    printf "%-18s" "$branch" && \
    timeout 20s ./ss -n 256 -f 6 -s 42 2>/dev/null | grep "FPS=" | tail -n +5 \
      | awk -F'= ' '{s+=$2; c++} END {printf "%.2f FPS\n", s/c}'
done
git checkout sequential >/dev/null 2>&1
```
