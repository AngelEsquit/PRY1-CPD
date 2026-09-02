# Fluid Screensaver — Navier-Stokes

Project #1 — Parallel and Distributed Computing, UVG.

A screensaver that simulates an incompressible 2D fluid by solving the
Navier-Stokes equations with Jos Stam's **Stable Fluids** method (1999).
Several pseudo-random colored ink sources inject color and momentum; the
fluid carries them around, generating swirls. The sources themselves move
under a mutual-gravity n-body system.

**This branch — `01-rb-tree`**: n-body gravity uses a Barnes-Hut quadtree
(`O(f log f)`) instead of brute-force pairwise gravity (`O(f²)`). Still
single-threaded, no OpenMP yet — this is a pure algorithmic step before
parallelization starts. ~10.51 FPS at the default config (n=256, f=6),
statistically the same as the `sequential` baseline at this resolution:
the fluid solver dominates the frame cost so much that the n-body
algorithm doesn't move the needle yet. See `docs/Roadmap.md` for the full
benchmark table.

---

## Branch model

This project is built in incremental steps, one branch per step, each
measured before merging into `master`:

| Branch | What it is |
|---|---|
| `sequential` | Frozen baseline: 1 thread, brute-force n-body, row-by-row Gauss-Seidel. Never touched again. |
| `master` | Integration branch. Always the best working version so far. |
| `01-rb-tree`, `02-omp-loops`, ... | One step each, branched from `master`, measured, then merged back. |
| `parallel-omp` | Frozen external reference (pre-roadmap OpenMP version, own legacy CLI). Never merged in. |

Full plan, benchmark methodology and results: `docs/Roadmap.md`.

---

## Requirements

- A C compiler with C11 support (`gcc` or `clang`)
- SDL2 (development library)
- `make`, `pkg-config`

### Installing SDL2

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

**Windows** — using WSL2 with Ubuntu and following the Debian instructions is
recommended. With MSYS2: `pacman -S mingw-w64-x86_64-SDL2 mingw-w64-x86_64-gcc`.

Verify the install:
```bash
pkg-config --modversion sdl2
```

---

## Building

```bash
make
```

Produces the `ss` executable (short for "screensaver"). To clean: `make clean`.

---

## Usage

```bash
./ss [options]
```

| Option | Description | Range | Default |
|--------|-------------|-------|---------|
| `-n <N>` | Grid resolution (N × N cells) | 16–1024 | 256 |
| `-f <F>` | Number of ink sources | 1–256 | 6 |
| `-W <width>` | Window width in pixels | ≥ 640 | 1920 |
| `-H <height>` | Window height in pixels | ≥ 480 | 1080 |
| `-s <seed>` | Pseudo-random seed (for reproducibility) | ≥ 0 | system clock |
| `-v <visc>` | Kinematic viscosity of the fluid | 0–1 | 0.0 |
| `-d <diff>` | Ink diffusion coefficient | 0–1 | 0.0001 |
| `-p`, `-F` | Fullscreen (exact size of the current screen) | — | off |
| `-h` | Shows help | — | — |

### Examples

```bash
./ss                          # default configuration
./ss -n 192 -f 12             # more resolution and more sources
./ss -n 256 -f 8 -s 42        # reproducible with a fixed seed
./ss -f 10                    # 10 sources, moving under mutual gravity
./ss -p                       # fullscreen at the current resolution
```

### Controls

| Key | Action |
|-----|--------|
| `ESC` / `Q` | Quit |
| `R` | Reset the simulation with new sources |

FPS are shown in the window title and printed to the console every 0.5 s.

---

## Algorithm structure

Each frame runs one time step made up of four operators:

1. **`add_source`** — applies external forces and ink injection.
2. **`diffuse`** — implicit viscous diffusion (sparse linear system).
3. **`advect`** — semi-Lagrangian transport with bilinear interpolation.
4. **`project`** — Hodge projection enforcing `div(u) = 0`
   (incompressibility); requires solving a Poisson equation.

Steps 2 and 4 solve linear systems with **Gauss-Seidel relaxation**
(`solve_linear`, `GAUSS_SEIDEL_ITERS` = 20 iterations).

Ink is carried in **three independent fields (R, G, B)** so colors can mix
on screen, which triples the cost of the density step.

**Trigonometric element:** each source's jet direction rotates over time as
`(cos(phase), sin(phase))` with its own angular velocity, generating spiral
vortices. Boundary conditions also reflect the normal velocity component,
making the fluid bounce off the walls.

---

## Parallelization in this branch

Still fully single-threaded — no OpenMP directives anywhere in this
checkout. The only change from `sequential` is the n-body algorithm
(Barnes-Hut quadtree instead of brute force), which is a pure
data-structure change, not a parallelization step.

`solve_linear` still uses row-by-row Gauss-Seidel, which reads values
**already updated within the same iteration**. That data dependency
between neighboring cells is what blocks a naive `#pragma omp parallel
for` here — see `docs/Roadmap.md` (Paso 3) for how later branches solve
it with a red-black checkerboard split.

---

## Verification performed

- Compiles without warnings under `-Wall -Wextra -O2 -std=c11`.
- Argument validation tested with non-numeric input, out-of-range values,
  missing values, and unknown options.
- `valgrind --leak-check=full`: 0 bytes lost, 0 errors.
