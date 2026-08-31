# Project development: structure and workflow

## 1. Project goal

This repository implements a 2D fluid screensaver based on the
Navier-Stokes equations for incompressible fluids, using Jos Stam's
"Stable Fluids" approach.

## 2. Current state of the repository

The project is a standalone, functional sequential implementation that
compiles cleanly. Parallelization is tracked separately via git branches
rather than a reserved subfolder: `sequential` is the baseline developed
here, and `parallel-omp` holds the OpenMP version. See
`docs/OPTIMIZATIONS.md` for the branch map, the step-by-step
parallelization plan, and benchmark results as they get filled in.

---

## 3. Project structure

```text
PRY1-CPD/
├── Makefile
├── README.md
├── ss
├── docs/
│   ├── Dev.md
│   └── OPTIMIZATIONS.md
├── include/
│   ├── background.h
│   ├── common.h
│   ├── config.h
│   ├── fields.h
│   ├── nbody.h
│   ├── render.h
│   ├── solver.h
│   ├── sources.h
│   └── utils.h
└── src/
    ├── background.c
    ├── config.c
    ├── fields.c
    ├── main.c
    ├── nbody.c
    ├── render.c
    ├── solver.c
    ├── sources.c
    └── utils.c
```

### Organization convention

- No `.c` files at the repo root.
- Each module lives directly under `src/`, one file per responsibility.
- Shared, reusable headers live in `include/`.
- Documentation lives in `docs/`.
- The built binary (`ss`) lands at the repo root, per the Makefile.

---

## 4. Functional description

The program simulates an incompressible 2D fluid flow over a regular grid.
The system is solved with a time-stepping approach:

1. Ink and force injection from random sources.
2. Diffusion of ink and velocity.
3. Advection of the field via semi-Lagrangian transport.
4. Hodge projection to keep the fluid incompressible.

The renderer visualizes the ink field as a texture shown in an SDL2 window.

Key variables:

- `grid_n`: grid resolution
- `sources`: number of ink injection points
- `dt`: time step
- `viscosity`: controls fluid damping
- `diffusion`: controls color spread

---

## 5. Responsibility of each module

### include/common.h
Project-wide global definitions:

- physical and validation constants
- the `IX()` macro for grid indexing
- boundary codes
- system default values

### include/config.h
Defines the program's config struct and the interface for parsing
command-line arguments:

- -n
- -f
- -W
- -H
- -s
- -v
- -d
- -p / -F
- -h

### include/fields.h
Defines the main struct holding the fluid's state:

- x/y velocities
- RGB ink fields
- scratch buffers for the numerical steps
- grid metadata

### include/sources.h
Defines the ink source model and how it's injected:

- position
- force
- angular direction
- ink flow rate

### include/nbody.h
Defines the n-body system that drives source movement (mutual gravity,
elastic bouncing off the grid edges).

### include/solver.h
Exposes the numerical core of the algorithm:

- diffusion
- advection
- projection
- buffer swapping
- Gauss-Seidel linear solve

### include/render.h
Defines the render API for filling the SDL texture and converting the ink
field's data into a visible image.

### include/background.h
Defines the twinkling starfield drawn behind the ink.

### include/utils.h
General-purpose utilities:

- value clamping
- random number generation

---

## 6. Current source files

### src/main.c
The program's entry point. Here it:

- parses arguments
- allocates field memory
- creates the SDL window
- runs the main loop
- updates the simulation per frame
- renders the texture
- measures FPS and updates the window title

### src/config.c
Implements validation and reading of the program's options. Keeps parsing
logic separate and isolated from the rest of the main flow.

### src/fields.c
Manages allocation, freeing and clearing of the fluid fields. This is the
system's memory layer.

### src/sources.c
Generates the ink sources and injects them during the simulation. Controls
the time variation and spin of the jet direction here.

### src/nbody.c
Drives the sources' movement as a mutual-gravity n-body system.

### src/solver.c
Implements the core numerics: diffusion, advection, boundary conditions,
projection and linear relaxation. This is the most important block for the
parallel version (see the OpenMP note in the file itself).

### src/render.c
Connects the numerical state to SDL to draw the simulation on screen.

### src/background.c
Renders the twinkling starfield behind the ink.

### src/utils.c
Centralizes helper functions that belong to neither the physics core nor
the visual layer.

---

## 7. Parallel version workflow

Parallelization work happens on its own git branch(es) built on top of
`sequential`, not in a separate source subfolder — see
`docs/OPTIMIZATIONS.md` for the branch map and the step-by-step
incremental plan (one optimization per branch, same benchmark suite rerun
after each step).

Recommended approach:

1. Keep the module API stable on `sequential`.
2. Implement parallel changes incrementally, one concern per branch,
   building on the previous step.
3. Benchmark under the same configuration after each step.
4. Document behavior differences and results in `docs/OPTIMIZATIONS.md`.

---

## 8. Building and running

From the project root:

```bash
make
./ss -n 128 -f 6
```

To clean:

```bash
make clean
```

The build uses SDL2 and gcc with C11 support.

---

## 9. Recommendations for the parallel phase

See `docs/OPTIMIZATIONS.md` for the concrete roadmap. In short:

1. Analyze the hot-spot kernels in `solver.c`.
2. Identify directly parallelizable elements: advection, injection,
   dissipation, the projection loops outside the linear solver.
3. Handle Gauss-Seidel relaxation carefully, since it introduces data
   dependencies between neighboring cells (see the note in `solver.c`).
4. Keep the same input/output API for a fair sequential-vs-parallel
   comparison.
5. Measure FPS and speedup with a reproducible benchmark suite.

---

## 10. Summary

The current structure meets the order and maintainability expected for a
parallelization phase:

- code organized by module, one file per responsibility
- centralized includes
- technical documentation in `docs/`
- a clean root with no `.c` files
- a functional, reproducible Makefile
- parallel work tracked via git branches, documented in
  `docs/OPTIMIZATIONS.md`

This leaves the project ready to evolve toward the parallel version without
losing traceability or clarity.
