# Robertson VODE Example

This repository provides a standalone build of the AMReX Microphysics Robertson stiff ODE example. It links against the `Microphysics` integration libraries (VODE, network definitions, interfaces) to evolve the classic three-species Robertson system on the CPU.

## Prerequisites

- CMake 3.20+ and a C++17 compiler
- Git and an internet connection (CMake fetches the pinned AMReX snapshot)
- A POSIX-like shell (commands are shown for macOS/Linux)

## Configure

The project uses out-of-source builds. Configure into `build/` (created automatically) and pin the default AMReX `development` tag:

```bash
cmake -S . -B build -DAMREX_TAG=development
```

You may set `-DAMREX_TAG=<tag_or_sha>` to test a specific AMReX revision. The configuration step downloads AMReX and wires in the Microphysics `interfaces`, `integration`, and `util` directories from the parent tree.

## Build

Compile the example target with:

```bash
cmake --build build --target robertson_vode -j
```

Adjust `-j` to match your core count. All generated artifacts stay under `build/`.

## Run

Execute the binary and capture its tabulated solution:

```bash
./build/robertson_vode > robertson.log
```

Inspect `robertson.log` to confirm the stiff solution curve (notably the small `y2` peak near `4e2`). When you update tolerances or chemistry, keep a reference log to compare future runs with `diff`.

## Plot

You can visualize the abundances with the helper script `plot_robertson.py` (requires `matplotlib`, e.g., `python3 -m pip install matplotlib`):

```bash
./build/robertson_vode > robertson.log   # if you have not already
python3 plot_robertson.py                # reads robertson.log by default
```

Pass an explicit log path (e.g., `python3 plot_robertson.py custom.log`) or save the figure instead of showing an interactive window with `-o output.png`.

## Repository Layout

- `src/`: user-editable C++ sources (currently `main.cpp`)
- `include/`: headers for integrator hooks, network definitions, and runtime parameters
- `build/`: generated files (not version-controlled)

Follow the AMReX/Microphysics coding style—four-space indentation, `CamelCase` types, `snake_case` functions/variables, and `.H` header extensions. Keep GPU annotations unchanged if you add device kernels.
