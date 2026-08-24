# Screensaver de Fluidos — Navier-Stokes (versión secuencial)

Proyecto #1 — Computación Paralela y Distribuida, UVG.

Screensaver que simula un fluido incompresible en 2D resolviendo las ecuaciones
de Navier-Stokes con el método **Stable Fluids** de Jos Stam (1999). Varias
fuentes de tinta de colores pseudoaleatorios inyectan color y cantidad de
movimiento; el fluido las transporta generando remolinos.

Esta es la **versión secuencial**, que sirve como base de comparación para medir
speedup y eficiencia de la versión paralela con OpenMP.

---

## Requisitos

- Compilador C con soporte C11 (`gcc` o `clang`)
- SDL2 (biblioteca de desarrollo)
- `make`, `pkg-config`

### Instalación de SDL2

**Ubuntu / Debian / WSL**
```bash
sudo apt-get update
sudo apt-get install build-essential pkg-config libsdl2-dev
```

**Fedora**
```bash
sudo dnf install gcc make pkgconf-pkg-config SDL2-devel
```

**macOS (Homebrew)**
```bash
brew install sdl2 pkg-config
```

**Windows** — se recomienda usar WSL2 con Ubuntu y seguir las instrucciones de
Debian. Si se usa MSYS2: `pacman -S mingw-w64-x86_64-SDL2 mingw-w64-x86_64-gcc`.

Verificar la instalación:
```bash
pkg-config --modversion sdl2
```

---

## Compilación

```bash
make
```

Genera el ejecutable `screensaver_seq`. Para limpiar: `make clean`.

---

## Uso

```bash
./screensaver_seq [opciones]
```

| Opción | Descripción | Rango | Defecto |
|--------|-------------|-------|---------|
| `-n <N>` | Resolución de la malla (N × N celdas) | 16–1024 | 128 |
| `-f <F>` | Cantidad de fuentes de tinta | 1–256 | 6 |
| `-W <ancho>` | Ancho de la ventana en píxeles | ≥ 640 | 800 |
| `-H <alto>` | Alto de la ventana en píxeles | ≥ 480 | 600 |
| `-s <semilla>` | Semilla pseudoaleatoria (para reproducibilidad) | ≥ 0 | reloj |
| `-t <dt>` | Paso de tiempo de la simulación | 0.001–1.0 | 0.10 |
| `-v <visc>` | Viscosidad cinemática del fluido | 0–1 | 0.0 |
| `-d <diff>` | Coeficiente de difusión de la tinta | 0–1 | 0.0 |
| `-b` | Activa el sistema de n-cuerpos (las fuentes se mueven por gravedad) | — | desactivado |
| `-h` | Muestra la ayuda | — | — |

### Ejemplos

```bash
./screensaver_seq                          # configuración por defecto
./screensaver_seq -n 192 -f 12             # más resolución y más fuentes
./screensaver_seq -n 256 -f 8 -s 42        # reproducible con semilla fija
./screensaver_seq -n 128 -v 0.0001 -d 0.00001   # fluido más viscoso y difuso
./screensaver_seq -b -f 10                 # fuentes en movimiento (n-cuerpos)
```

### Controles

| Tecla | Acción |
|-------|--------|
| `ESC` / `Q` | Salir |
| `R` | Reiniciar la simulación con fuentes nuevas |

Los FPS se muestran en el título de la ventana y se imprimen en consola cada
0.5 s.

---

## Estructura del algoritmo

Cada frame ejecuta un paso de tiempo compuesto por cuatro operadores:

1. **`agregar_fuente`** — aplica fuerzas externas y la inyección de tinta.
2. **`difundir`** — difusión viscosa implícita (sistema lineal disperso).
3. **`advectar`** — transporte semi-Lagrangiano con interpolación bilineal.
4. **`proyectar`** — proyección de Hodge que impone `div(u) = 0`
   (incompresibilidad); requiere resolver una ecuación de Poisson.

Los pasos 2 y 4 resuelven sistemas lineales con **relajación de Gauss-Seidel**
(`resolver_lineal`, 20 iteraciones).

La tinta se lleva en **tres campos independientes (R, G, B)** para poder mezclar
colores en pantalla, lo que triplica el costo del paso de densidad.

**Elemento trigonométrico:** la dirección del chorro de cada fuente rota en el
tiempo según `(cos(fase), sin(fase))` con velocidad angular propia, generando
vórtices en espiral. Además, las condiciones de frontera reflejan la componente
normal de la velocidad, haciendo que el fluido rebote en los muros.

---

## Nota para la fase paralela

`resolver_lineal` usa Gauss-Seidel, que lee valores **ya actualizados en la misma
iteración**. Esto crea una dependencia de datos entre celdas vecinas y hace que
el ciclo **no sea directamente paralelizable** con `#pragma omp parallel for`
(produciría condiciones de carrera y resultados no deterministas).

Opciones para la versión OpenMP:

- **Jacobi** — leer de un buffer separado elimina las dependencias por completo;
  es trivialmente paralelizable pero converge más lento por iteración (puede
  requerir aumentar `ITER_GAUSS_SEIDEL`).
- **Red-black Gauss-Seidel** — dividir la malla en celdas "rojas" y "negras"
  como un tablero de ajedrez; cada mitad se actualiza en paralelo sin
  dependencias internas, conservando la velocidad de convergencia de
  Gauss-Seidel.

Los operadores `advectar`, `agregar_fuente`, `disipar_tinta` y los ciclos de
`proyectar` (fuera del solver) **sí** son directamente paralelizables, ya que
cada celda escribe únicamente en su propia posición.

---

## Mediciones de referencia (versión secuencial)

Medido con `-f 8 -s 7`, promedio de las últimas 3 lecturas de FPS:

| N (malla) | Celdas | FPS medio |
|-----------|--------|-----------|
| 64 | 4 096 | 154.2 |
| 128 | 16 384 | 51.0 |
| 192 | 36 864 | 23.9 |
| 256 | 65 536 | 13.7 |

> Estos valores dependen del hardware; cada equipo debe repetir las mediciones
> en su propia máquina para la bitácora de pruebas.

**Observación relevante para el informe:** a partir de N ≈ 192 la versión
secuencial cae por debajo del umbral de 30 FPS exigido por el enunciado. Esto
justifica directamente la necesidad de la versión paralela y da un punto de
comparación muy claro para la sección de speedup.

---

## Verificación realizada

- Compila sin warnings con `-Wall -Wextra -O2 -std=c11`.
- Validación de argumentos probada con entradas no numéricas, fuera de rango,
  faltantes y opciones desconocidas.
- `valgrind --leak-check=full`: 0 bytes perdidos, 0 errores.