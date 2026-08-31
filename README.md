# Fluid Screensaver — Navier-Stokes (sequential version)

Project #1 — Parallel and Distributed Computing, UVG.

A screensaver that simulates an incompressible 2D fluid by solving the
Navier-Stokes equations with Jos Stam's **Stable Fluids** method (1999).
Several pseudo-random colored ink sources inject color and momentum; the
fluid carries them around, generating swirls. The sources themselves move
under a mutual-gravity n-body system.

This is the **sequential version**, which serves as the comparison baseline
for measuring speedup and efficiency of the OpenMP parallel version (see
`docs/OPTIMIZATIONS.md` for the parallelization plan and benchmark results).

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

Produces the `Screensaver` executable. To clean: `make clean`.

---

## Usage

```bash
./Screensaver [options]
```

| Option | Description | Range | Default |
|--------|-------------|-------|---------|
| `-n <N>` | Grid resolution (N × N cells) | 16–1024 | 256 |
| `-f <F>` | Number of ink sources | 1–256 | 6 |
| `-W <width>` | Window width in pixels | ≥ 640 | 1920 |
| `-H <height>` | Window height in pixels | ≥ 480 | 1080 |
| `-s <seed>` | Pseudo-random seed (for reproducibility) | ≥ 0 | system clock |
| `-t <dt>` | Simulation time step | 0.001–1.0 | 0.07 |
| `-v <visc>` | Kinematic viscosity of the fluid | 0–1 | 0.0 |
| `-d <diff>` | Ink diffusion coefficient | 0–1 | 0.0001 |
| `-b` | Toggles the n-body system (sources move under mutual gravity) | — | on |
| `-p`, `-F` | Fullscreen (exact size of the current screen) | — | off |
| `-h` | Shows help | — | — |

### Examples

```bash
./Screensaver                          # default configuration
./Screensaver -n 192 -f 12             # more resolution and more sources
./Screensaver -n 256 -f 8 -s 42        # reproducible with a fixed seed
./Screensaver -n 128 -v 0.0001 -d 0.00001   # more viscous, more diffuse fluid
./Screensaver -b -f 10                 # disable n-body motion, 10 sources
./Screensaver -p                       # fullscreen at the current resolution
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

## Note on the parallel phase

`solve_linear` uses Gauss-Seidel, which reads values **already updated
within the same iteration**. This creates a data dependency between
neighboring cells and makes the loop **not directly parallelizable** with
`#pragma omp parallel for` (it would produce race conditions and
non-deterministic results).

Options for the OpenMP version:

- **Jacobi** — reading from a separate buffer removes the dependencies
  entirely; trivially parallelizable but converges slower per iteration
  (may require raising `GAUSS_SEIDEL_ITERS`).
- **Red-black Gauss-Seidel** — split the grid into "red" and "black" cells
  like a checkerboard; each half updates in parallel with no internal
  dependencies, keeping Gauss-Seidel's convergence speed.

The `advect`, `add_source`, `dissipate_ink` operators and the loops inside
`project` (outside the solver) **are** directly parallelizable, since each
cell only writes to its own position.

See `docs/OPTIMIZATIONS.md` for the actual step-by-step parallelization plan
and benchmark results as they're filled in.

---

## Verification performed

- Compiles without warnings under `-Wall -Wextra -O2 -std=c11`.
- Argument validation tested with non-numeric input, out-of-range values,
  missing values, and unknown options.
- `valgrind --leak-check=full`: 0 bytes lost, 0 errors.
