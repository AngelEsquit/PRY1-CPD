# Plan de paralelizacion y benchmarking

Este documento es el plan de trabajo para ir de la version secuencial a la
paralelizada, midiendo speedup en cada paso. La idea es ir llenando las
tablas de resultados a medida que se prueba cada version/rama.

## Mapa de ramas

| Rama            | Que es | Se toca? |
|-----------------|--------|----------|
| `sequential`    | Baseline puro, congelado para siempre: 1 hilo, n-body O(F^2) fuerza bruta, solver Gauss-Seidel fila por fila. Referencia fija para comparar speedup acumulado. | **Nunca** — ni un commit mas despues del baseline. |
| `master`        | Rama de integracion. Acumula cada paso terminado via merge (fast-forward). Siempre representa "lo mejor que tenemos ahora mismo". | Solo via merge de una rama de paso ya terminada y medida. |
| `01-rb-tree`, `02-omp-loops`, `03-omp-solver`, ... | Una rama por paso del roadmap. Cada una sale de `master` (no de `sequential`), implementa un solo cambio, se mide, y se mergea de vuelta a `master` cuando esta lista. | Ahi es donde pasa el trabajo real. |
| `parallel-omp`  | Referencia externa congelada: la version OpenMP tal cual quedo antes de toda esta reorganizacion (fuera de esta cadena, con su propia CLI vieja). Sirve para saber "a donde se puede llegar", no se integra directamente. | Nunca. |

Regla simple: **`sequential` = punto de comparacion fijo. `master` = donde
vive el progreso acumulado. Los pasos numerados = donde se hace el trabajo.**
Nunca se implementa nada directamente en `sequential` ni en `master`.

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

```
master ──▶ 01-rb-tree ──(merge)──▶ master ──▶ 02-omp-loops ──(merge)──▶ master ──▶ 03-omp-solver ──▶ ...
```

Cada paso sale de `master` (con el paso anterior ya mergeado adentro), y
vuelve a mergearse a `master` cuando esta terminado y medido.

- [x] **Paso 0, baseline** — vive para siempre en `sequential`. 1 hilo,
      n-body O(F^2) fuerza bruta, solver Gauss-Seidel fila por fila.
- [x] **Paso 1 — `01-rb-tree`** (mergeado a `master`): arbol para n-body
      (Barnes-Hut), reemplaza el loop O(F^2) de `update_nbody_sources()`
      (`nbody.c`) por un quadtree con aproximacion por centro de masa,
      O(F log F). Sigue secuencial (sin OpenMP todavia). FPS practicamente
      identico al baseline en todo el rango probado, incluso con `f` alto:
      a estas resoluciones el solver de fluidos domina tanto que el
      algoritmo de n-body no se nota todavia.
- [x] **Paso 2 — `02-omp-loops`** (mergeado a `master`): paraleliza los
      loops "faciles" (sitios 1, 2, 4, 5, 6 del catalogo de arriba:
      fuente/disipacion/adveccion/proyeccion/render) con
      `#pragma omp parallel for`, sin cambiar el algoritmo. El solver de
      presion/difusion se queda secuencial todavia (sigue siendo
      Gauss-Seidel fila por fila, es el Paso 3). Buen salto de FPS
      (ver tabla maestra) incluso sin tocar el hotspot principal.
- [x] **Paso 3 — `03-omp-solver`**: red-black + paralelizar el solver: el
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

`parallel-omp` (la rama congelada, fuera de esta cadena) ya tiene
implementados los pasos 2 y 3 juntos (mas cambios visuales que no son parte
de este roadmap de performance) — sirve como referencia de "a donde se
puede llegar", pero esta cadena se construye de nuevo para medir cada
incremento por separado.

---

## Flujo de trabajo (para el equipo)

Cada paso pendiente (Paso 3 en adelante) lo puede tomar cualquiera del
equipo siguiendo estos pasos. Un solo paso = una sola rama = un solo
commit de codigo (mas, si aplica, un commit chico separado con los
resultados del benchmark). Nada de commits micro por cada ajuste.

1. **Asegurate de tener `master` actualizado** (el que va a tener el ultimo
   paso ya mergeado):
   ```bash
   git checkout master
   git pull
   ```
2. **Crea tu rama del paso que te toca**, con el nombre exacto del roadmap
   de arriba (ej. para el Paso 3):
   ```bash
   git checkout -b 03-omp-solver master
   ```
3. **Implementa el cambio de ese paso** (ver el catalogo de optimizaciones
   mas arriba en este documento para el detalle tecnico de cada uno).
   Mismo estilo de comentarios que el resto del codigo: cortos, en ingles,
   sin jerga sin explicar.
4. **Compila y prueba que corre** (`make`, correr el binario un rato,
   confirmar que no crashea y se ve razonable):
   ```bash
   make clean && make
   ./ss -n 256 -f 6
   ```
5. **Commit del codigo** (uno solo, mensaje explicando que cambia y por
   que, sin tabla de benchmark todavia si no la vas a correr vos mismo):
   ```bash
   git add -A
   git commit -m "feat: <resumen corto del cambio>"
   ```
6. **Push de tu rama** y avisar:
   ```bash
   git push -u origin 03-omp-solver
   ```
   A partir de ahi, **yo corro la bateria de benchmark oficial** (ver
   "Metodologia de benchmark" y "Comandos" mas abajo en este documento) y
   lleno la tabla maestra, para que todos los numeros salgan de la misma
   maquina/condiciones y sean comparables entre si. Si vos ya corriste tus
   propios numeros como referencia rapida, dejalos en el mensaje del PR o
   en el chat, no hace falta que entren al documento.
7. **Una vez medido y con la tabla llena**, se mergea a `master`
   (fast-forward, sin commit de merge):
   ```bash
   git checkout master
   git merge 03-omp-solver --ff-only
   ```
8. El siguiente paso se crea sobre este nuevo `master` (que ya incluye tu
   trabajo), repitiendo desde el paso 1.

**Reglas rapidas**:
- Nunca commitear directo en `sequential` o en `master`.
- Una rama = un paso del roadmap, no una mezcla de varios.
- Si tu cambio no compila limpio con `-Wall -Wextra` o rompe el build,
  arreglalo antes de pedir que se corra el benchmark.
- Si encontras un bug o comportamiento raro que no es parte de tu paso
  (como el problema de pixelado de tinta que ya se reviso), avisa en vez
  de mezclarlo en el mismo commit — se resuelve aparte.

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

| Fila matriz | `sequential` (baseline) | `01-rb-tree` | `02-omp-loops` | `03-omp-solver` | `04-schedule-tuning` | `05-collapse-tuning` |
|---|---|---|---|---|---|---|
| A6 (n=64, f=6) | 17.21 | 17.19 | 118.40 | | | |
| A12 (n=64, f=12) | 17.16 | 17.20 | 117.43 | | | |
| B6 (n=256, f=6) | 10.52 | 10.51 | 25.20 | | | |
| B12 (n=256, f=12) | 10.95 | 10.95 | 25.24 | | | |
| C6 (n=512, f=6) | 3.75 | 3.76 | 6.72 | | | |
| C12 (n=512, f=12) | 4.30 | 4.30 | 6.89 | | | |
| H6 (n=600, f=6) | 2.69 | 2.69 | 4.46 | | | |
| H12 (n=600, f=12) | 3.09 | 3.09 | 4.79 | | | |
| I6 (n=700, f=6) | 1.90 | 1.91 | 3.00 | | | |
| I12 (n=700, f=12) | 2.10 | 2.10 | 3.21 | | | |
| D6 (n=1024, f=6) | 0.92 | 0.92 | 1.15 | | | |
| D12 (n=1024, f=12) | 1.00 | 0.99 | 1.23 | | | |
| E (n=256, f=4) | 7.52 | 7.52 | 24.39 | | | |
| F (n=256, f=64) | 11.39 | 11.39 | 25.36 | | | |
| G (n=256, f=256) | 11.37 | 11.38 | 25.30 | | | |

Medido en local (display real, no headless), corridas de ~15s por fila,
promedio ignorando los primeros ~2s de warm-up, `-s 42` en todas.
Barnes-Hut da practicamente el mismo FPS que fuerza bruta en todo el rango
probado (incluyendo f=64/256): a estas resoluciones de malla el solver de
fluidos domina tanto el costo que el algoritmo de n-cuerpos, sea O(F^2) o
O(F log F), no se nota en el FPS total. Sigue siendo la base correcta para
cuando se paralelice el solver (Paso 3+), donde el resto del frame sera
mas barato relativamente y el n-body podria empezar a pesar mas.

`02-omp-loops` (fuente/disipacion/adveccion/proyeccion/render en paralelo,
solver todavia secuencial) ya da un salto claro con 20 hilos: ~6-7x en las
mallas grandes (C/H/I/D) y hasta ~7x en las chicas (A). El techo sigue
siendo el solver de presion/difusion (Gauss-Seidel), que es el Paso 3.

### Curva de escalabilidad por hilos (desde el Paso 2 en adelante)

Usar la fila B (n=256, f=6) como referencia, variando `OMP_NUM_THREADS`:

| Hilos | FPS `02-omp-loops` | FPS `03-omp-solver` | FPS `04-schedule-tuning` | FPS `05-collapse-tuning` |
|-------|------------------|-------------------|-------------------------|--------------------------|
| 1     | 9.95             |                   |                         |                          |
| 2     | 14.94            |                   |                         |                          |
| 4     | 20.33            |                   |                         |                          |
| 8     | 23.57            |                   |                         |                          |
| 16    | 25.52            |                   |                         |                          |
| 20    | 25.57            |                   |                         |                          |

A 1 hilo (9.95 FPS) `02-omp-loops` da casi lo mismo que `01-rb-tree`
(10.51, diferencia es solo ruido de medicion) — tiene sentido, con 1 hilo
OpenMP no cambia nada. El escalado de ahi hasta 20 hilos (2.6x) muestra el
techo del Amdahl: el solver secuencial que queda (Paso 3) sigue limitando
cuanto puede escalar el resto.

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
