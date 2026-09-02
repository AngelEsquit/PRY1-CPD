# Fluid Screensaver — Navier-Stokes

Project #1 — Parallel and Distributed Computing, UVG.

A screensaver that simulates an incompressible 2D fluid by solving the
Navier-Stokes equations with Jos Stam's **Stable Fluids** method (1999).
Several pseudo-random colored ink sources inject color and momentum; the
fluid carries them around, generating swirls. The sources themselves move
under a mutual-gravity n-body system.

**This branch — `05-collapse-tuning`**: was meant to compare `collapse(2)`
on/off on the advect/project/render loops, especially for small grids
(fewer rows than threads). **Cancelled before reliable numbers came out**
— the code here is functionally identical to `03-omp-solver`
(`schedule(dynamic, 16)`, `collapse(2)` still in place, unchanged). Not
pursued further; see `docs/Roadmap.md` for what happened. The performance
work that stuck landed in `04-schedule-tuning` instead.

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
- `make`, `pkg-config`, OpenMP support in the compiler (standard in gcc/clang)

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

Control the thread count without recompiling via the standard OpenMP
environment variable:
```bash
OMP_NUM_THREADS=4 ./ss -n 256 -f 6
```

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

Same as `03-omp-solver`: everything parallel (source injection,
dissipation, advection, projection, render, red-black solver), all with
`schedule(dynamic, 16)` and `collapse(2)` on the advect/project/render
loops. No changes were actually kept on this branch — the planned
collapse-on/off comparison never produced reliable numbers before being
cancelled.

---

## Verification performed

- Compiles without warnings under `-Wall -Wextra -O2 -std=c11 -fopenmp`.
- Argument validation tested with non-numeric input, out-of-range values,
  missing values, and unknown options.
- `valgrind --leak-check=full`: 0 bytes lost, 0 errors.
