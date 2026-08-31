# Plan de paralelizacion y benchmarking

Este documento es el plan de trabajo para ir de la version secuencial a la
paralelizada, midiendo speedup en cada paso. La idea es ir llenando las
tablas de resultados a medida que se prueba cada version/rama.

## Mapa de ramas

| Rama            | Binario        | Hilos | Solver de presion/difusion      | Fullscreen |
|-----------------|----------------|-------|----------------------------------|------------|
| `master`        | `Screensaver`  | 1     | Gauss-Seidel fila por fila       | si (`-p`)  |
| `sequential`    | `Screensaver`  | 1     | Gauss-Seidel fila por fila       | si (`-p`)  |
| `parallel-omp`  | `Screensaver`  | N     | Red-black (paralelizable)        | no         |

`sequential` = referencia correcta para medir speedup: mismos parametros,
mismo comportamiento visual, sin OpenMP y con el algoritmo de solver que de
verdad es secuencial (no el red-black, que solo existe para poder
paralelizarse — ver la nota en `solver.c`).

`parallel-omp` es la version con OpenMP tal cual quedo tras el refactor de
render/fuentes/fondo, con `-fopenmp` en el Makefile. **A partir de aqui se
deja como referencia congelada** (para ver el resultado final de a donde se
quiere llegar) y todo el trabajo de optimizacion incremental se hace sobre
`sequential`, creando una rama nueva por cada paso del checklist mas abajo.
Por eso `parallel-omp` se quedo con `-b` desactivado por default y sin
fullscreen — no se le sigue actualizando a la par de `sequential`. Al
comparar FPS entre ramas, pasar los flags explicitamente (`-b`) en vez de
confiar en el default, ya que default difiere entre ramas.

Cuando se cree una rama nueva por cada optimizacion incremental (p.ej.
`parallel-v1-omp-basico`, `parallel-v2-schedule-tuning`, etc.), agregar una
fila a la tabla de arriba y una seccion nueva en "Bitacora de optimizaciones"
mas abajo.

---

## Parametros (linea de comandos y constantes)

### Flags de `./Screensaver`

| Flag | Rango / default | Efecto |
|------|------------------|--------|
| `-n <N>` | `[16..1024]`, def. `256` | Resolucion de la malla (N x N celdas interiores). El costo del solver escala ~O(N^2) por iteracion; es el parametro mas importante para medir speedup. |
| `-f <F>` | `[1..256]`, def. `6` | Cantidad de fuentes de tinta. |
| `-W <ancho>` | min. `640`, def. `1920` | Ancho de ventana en px (no afecta el costo del solver, solo el de render). |
| `-H <alto>` | min. `480`, def. `1080` | Alto de ventana en px. |
| `-s <semilla>` | def. reloj del sistema | Semilla PRNG — fijarla para comparaciones reproducibles entre corridas. |
| `-t <dt>` | `[0.001..1.0]`, def. `0.07` | Paso de tiempo de la simulacion. |
| `-v <visc>` | `[0.0..1.0]`, def. `0.0` | Viscosidad cinematica. |
| `-d <diff>` | `[0.0..1.0]`, def. `0.0001` | Coeficiente de difusion de tinta (el "blur"/suavizado). |
| `-b` | activado por default en `sequential` (desactivado en `parallel-omp`, congelada) | Alterna el sistema de n-cuerpos (fuentes con gravedad mutua). |
| `-p`, `-F` | desactivado por default | Pantalla completa (solo en `sequential`/`master` por ahora). |
| `-h` | — | Ayuda. |

### Constantes (`include/comun.h`, `include/fuentes.h`, `src/secuencial/ncuerpos/ncuerpos.c`)

Estas no son ajustables por CLI, pero son relevantes si se quiere variar el
"blur"/comportamiento en experimentos mas finos (requiere recompilar):

| Constante | Valor | Que hace |
|---|---|---|
| `ITER_GAUSS_SEIDEL` | `20` | Iteraciones del solver lineal por llamada — sube linealmente el costo de `resolver_lineal()`, que es el hotspot principal. |
| `DISIPACION` | `0.995` | Retencion de tinta por frame. |
| `BRILLO_FACTOR` / `CONTRASTE_FACTOR` / `SATURACION_FACTOR` | `0.6` / `1.0` / `2.4` | Solo afectan el render (color), no el costo computacional. |
| `FUENTE_RADIO_BASE` / `FUENTE_SIGMA_BASE` | `2` / `1.1` | Tamano del vecindario de inyeccion gaussiana por fuente; escala con `-n` via `fuente_escala()`. |
| `NBODY_G` / `NBODY_SUAVIZADO` / `NBODY_VEL_MAX` / `NBODY_RESTITUCION` | `0.3` / `6.0` / `5.0` / `0.9` | Fisica del sistema de n-cuerpos (solo con `-b`); el costo es O(F^2) en cantidad de fuentes. |

**Parametros clave para el plan de benchmarking**: `-n` (resolucion de malla,
el que mas importa), `-f` (cantidad de fuentes), `-b` (activa costo O(F^2)
extra), y el numero de hilos OpenMP (`OMP_NUM_THREADS`, ver abajo — no es un
flag del programa, es una variable de entorno).

---

## Catalogo de optimizaciones (parallel-omp)

Cada fila es un sitio de paralelizacion real en el codigo. "Como
implementarlo" describe el patron para poder reproducirlo si se hace una
rama incremental que las va agregando una por una.

### 1. `agregar_fuente()` — `src/secuencial/fluidos/fuentes.c`

- **Que hace**: suma el termino fuente (velocidad/tinta inyectada) a cada
  celda, `O(N^2)`, sin dependencias entre celdas.
- **Como se paralelizo**: `#pragma omp parallel for schedule(dynamic, 16)`
  sobre el loop plano de `celdas_total`.
- **Por que funciona sin cambios de algoritmo**: cada celda se escribe una
  sola vez y no depende de otras celdas — paralelizable "tal cual".
- **Como implementarlo**: agregar el pragma justo antes del `for`. Nada mas
  cambia.

### 2. `disipar_tinta()` — `src/secuencial/fluidos/fuentes.c`

- Mismo patron que (1): `#pragma omp parallel for schedule(dynamic, 16)`
  sobre un loop plano e independiente por celda (`tinta *= DISIPACION`).

### 3. Solver de presion/difusion: Gauss-Seidel -> red-black — `src/secuencial/fluidos/solver.c`

- **Que hace**: resuelve el sistema lineal disperso de la difusion implicita
  y de la proyeccion de Hodge, iterativamente (`ITER_GAUSS_SEIDEL = 20`
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
  1. Separar el cuerpo del loop en una funcion `relajar_color(campo,
     campo_previo, a, inverso_c, paridad)` que recorre solo las celdas de
     una paridad (`inicio = ((1+j)%2 == paridad) ? 1 : 2`, luego `i += 2`).
  2. Envolver las `ITER_GAUSS_SEIDEL` iteraciones en **una sola** region
     `#pragma omp parallel` (abrir/cerrar el team de hilos tiene costo fijo;
     hacerlo una vez en vez de 2*20 veces importa en mallas chicas).
  3. Dentro, cada `relajar_color()` usa `#pragma omp for schedule(dynamic,
     16)` (barrera implicita automatica al final).
  4. Aplicar la condicion de frontera con `#pragma omp single` (un solo hilo
     la aplica mientras los demas esperan en la barrera implicita del
     `single`).
- **Nota de correctitud numerica**: red-black converge a un resultado
  practicamente identico a Gauss-Seidel secuencial (de hecho a veces
  converge un poco mas rapido), asi que el comportamiento visual no deberia
  notarse distinto entre `sequential` y `parallel-omp`.

### 4. `advectar()` — `src/secuencial/fluidos/solver.c`

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

### 5. Proyeccion de Hodge (divergencia + gradiente) — `src/secuencial/fluidos/solver.c`

- Dos loops mas con el mismo patron que (4): cada celda destino es
  independiente. `#pragma omp parallel for collapse(2) schedule(dynamic,
  16)` en cada uno.

### 6. `renderizar_tinta()` — `src/secuencial/render/render.c`

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

## Optimizaciones pendientes / ideas para ramas incrementales

Para ir documentando "de X a Y, speedup Z%" por rama, algunas ideas de pasos
intermedios entre `sequential` y `parallel-omp` completo (crear una rama por
cada uno, medir, anotar en la bitacora):

- [ ] **v1**: solo paralelizar los loops "vergonzosamente paralelos" (1, 2,
      4, 5, 6 de la lista de arriba), dejar el solver de Gauss-Seidel
      secuencial (bottleneck real, pero mide cuanto dan los loops faciles
      solos).
- [ ] **v2**: agregar red-black + paralelizacion del solver (el cambio que
      mas deberia importar, dado que es el hotspot).
- [ ] **v3**: variar `schedule()` — probar `static` vs `dynamic` vs
      `guided`, y distintos tamanos de chunk (4, 16, 64), medir cual da
      mejor FPS a cada resolucion.
- [ ] **v4**: probar sin `collapse(2)` (paralelizar solo el loop externo
      `j`) contra con `collapse(2)`, especialmente relevante en mallas
      chicas (`-n 16` a `-n 64`, menos filas que `nproc`).
- [ ] **v5** (opcional, mas alla de OpenMP for-loops): SIMD explicito
      (`#pragma omp simd`) en los loops internos del solver, o revisar si
      el compilador ya autovectoriza con `-O2`/`-O3` (comparar `-O2` vs
      `-O3` como variable de build, no solo de runtime).

---

## Plan de experimentacion / benchmarking

### Metodologia

1. **Metrica**: el programa ya imprime FPS por stdout cada ~0.5s
   (`printf("FPS= %.2f\n", ...)`). Para un numero estable, correr cada
   configuracion un tiempo fijo (p.ej. 30s) y promediar las lineas de FPS
   ignorando las primeras ~5 (warm-up / primeros frames con cache fria).
2. **Reproducibilidad**: fijar siempre la semilla (`-s 42`) para que todas
   las corridas simulen exactamente las mismas fuentes/trayectorias.
3. **Automatizar la corrida** (ejemplo, ajustar a lo que se necesite):
   ```bash
   timeout 30s ./Screensaver -n 256 -f 6 -s 42 2>/dev/null \
     | grep "FPS=" | tail -n +6 \
     | awk -F'= ' '{s+=$2; c++} END {printf "FPS promedio: %.2f\n", s/c}'
   ```
   (En una maquina sin display, correr con `xvfb-run` o similar si hace
   falta un X server virtual.)
4. **Hilos OpenMP**: controlar la cantidad de hilos con la variable de
   entorno `OMP_NUM_THREADS`, sin recompilar:
   ```bash
   OMP_NUM_THREADS=4 ./Screensaver -n 256 -f 6 -s 42
   ```
   Nucleos disponibles en esta maquina: `nproc` → **20**.
5. **Speedup**: `speedup = FPS_paralelo / FPS_secuencial` (a igualdad de
   parametros de simulacion). **Eficiencia**: `speedup / num_hilos`.

### Tabla 1 — Speedup vs. resolucion de malla (`-n`)

Parametros fijos: `-f 6 -s 42`, sin `-b`. Correr `sequential` una vez
(no depende de hilos) y `parallel-omp` con `OMP_NUM_THREADS` = numero de
nucleos.

| `-n` | FPS `sequential` | FPS `parallel-omp` | Speedup | Notas |
|------|-------------------|----------------------|---------|-------|
| 32   |                   |                      |         |       |
| 64   |                   |                      |         |       |
| 128  |                   |                      |         |       |
| 256  |                   |                      |         |       |
| 512  |                   |                      |         |       |
| 1024 |                   |                      |         |       |

*Hipotesis a confirmar*: el speedup deberia crecer con `-n` — mallas chicas
no tienen suficiente trabajo para amortizar el overhead de sincronizacion
de OpenMP (barreras, creacion del team de hilos), mallas grandes si.

### Tabla 2 — Speedup vs. numero de hilos (`OMP_NUM_THREADS`)

Parametros fijos: `-n 512 -f 6 -s 42` (una malla grande, para que el
resultado no este dominado por overhead de sincronizacion).

| Hilos | FPS | Speedup vs. 1 hilo | Eficiencia (speedup/hilos) |
|-------|-----|---------------------|------------------------------|
| 1     |     |                     |                              |
| 2     |     |                     |                              |
| 4     |     |                     |                              |
| 8     |     |                     |                              |
| 16    |     |                     |                              |
| 20    |     |                     |                              |

*Hipotesis a confirmar*: la eficiencia deberia caer al acercarse al numero
de nucleos fisicos (ley de Amdahl — la parte secuencial del programa, como
la inicializacion por frame y el `single` de la frontera, pone un techo).

### Tabla 3 — Speedup vs. cantidad de fuentes (`-f`)

Parametros fijos: `-n 256 -s 42`, con y sin `-b` (n-cuerpos, O(F^2)).

| `-f` | FPS `sequential` (sin `-b`) | FPS `parallel-omp` (sin `-b`) | FPS `sequential` (con `-b`) | FPS `parallel-omp` (con `-b`) |
|------|-------------------------------|----------------------------------|--------------------------------|-----------------------------------|
| 4    |                               |                                   |                                 |                                    |
| 16   |                               |                                   |                                 |                                    |
| 64   |                               |                                   |                                 |                                    |
| 256  |                               |                                   |                                 |                                    |

*Nota*: `actualizar_fuentes_nbody()` (el sistema de n-cuerpos) **no esta
paralelizado** en ninguna rama actual — es O(F^2) pero solo importa con `-f`
alto. Buen candidato a rama incremental si se quiere otro punto de datos
("optimizacion 7").

### Tabla 4 — Bitacora de optimizaciones incrementales

Ir agregando una fila por cada rama nueva creada a partir del checklist de
"Optimizaciones pendientes":

| Rama | Que cambia respecto a la anterior | FPS (n=256,f=6, N hilos) | Speedup acumulado | Speedup incremental |
|------|-------------------------------------|-----------------------------|----------------------|-------------------------|
| `sequential` | baseline | | 1.00x | — |
| `parallel-omp` (v1: solo loops triviales) | | | | |
| `parallel-omp` (v2: + solver red-black) | | | | |
| ... | | | | |

---

## Como correr una comparacion rapida entre ramas

```bash
git checkout sequential  && make clean && make && \
  timeout 20s ./Screensaver -n 256 -f 6 -s 42 2>/dev/null | grep "FPS=" | tail -n +5

git checkout parallel-omp && make clean && make && \
  timeout 20s ./Screensaver -n 256 -f 6 -s 42 2>/dev/null | grep "FPS=" | tail -n +5
```
