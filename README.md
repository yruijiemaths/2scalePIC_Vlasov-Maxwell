# Two-Scale PIC Methods for the 1D--2V Vlasov--Maxwell System

This repository contains MPI-parallel Fortran implementations of uniformly
accurate two-scale particle-in-cell (PIC) methods for the
one-dimensional-in-space, two-dimensional-in-velocity (1D--2V)
Vlasov--Maxwell system under a strong external magnetic field.

The code provides four time integrators:

| CMake method | Description | Spatial discretization |
|---|---|---|
| `TSI-Sim` | symmetric implicit two-scale exponential integrator | collocated grid |
| `TSI-Sex` | symmetric explicit two-scale exponential integrator | collocated grid |
| `TSI-CPex` | charge-preserving explicit two-scale exponential integrator | collocated grid |
| `TSFV-CPimex` | charge-preserving two-scale finite-volume IMEX integrator | staggered grid |

The fast two-scale variable is discretized spectrally using FFTW. Particle
deposition and field interpolation use sixth-order B-spline shape functions.
The charge-preserving finite-volume method uses a compatible staggered-grid
current deposition.

## Requirements

- CMake 3.18 or newer;
- a Fortran compiler;
- MPI with Fortran bindings;
- FFTW3, including the `fftw3.f` header and FFTW3 library.

The project is intended for Linux/HPC environments. GNU Fortran and Intel
Fortran toolchains can be used with compatible MPI and FFTW builds.

## Directory structure

```text
PICfortran/
|-- config/
|   |-- pic.nml                    default runtime configuration
|   |-- example-two-stream.nml     small two-stream example
|   |-- example-weibel.nml         small Weibel example
|   `-- ...                        additional experiment configurations
|
|-- src/
|   |-- common/
|   |   |-- pic_config.F90         namelist loading and validation
|   |   `-- pic_output.F90         shared snapshots and diagnostics
|   |
|   |-- grids/
|   |   |-- standard/              collocated spatial discretization
|   |   |   |-- mesh.F90           mesh and field data types
|   |   |   |-- particles.F90      particle allocation and initialization
|   |   |   |-- interpolations.F90 B-spline interpolation and deposition
|   |   |   |-- poisson.F90        Poisson solver and spatial derivatives
|   |   |   `-- grid_output_adapter.F90
|   |   |                            grid-specific diagnostic operations
|   |   |
|   |   `-- half-grid/             staggered finite-volume discretization
|   |       |-- mesh.F90
|   |       |-- particles.F90
|   |       |-- interpolations.F90
|   |       |-- poisson.F90
|   |       `-- grid_output_adapter.F90
|   |
|   `-- methods/
|       |-- TSI-Sim/main.F90       symmetric implicit method
|       |-- TSI-Sex/main.F90       symmetric explicit method
|       |-- TSI-CPex/main.F90      charge-preserving explicit method
|       `-- TSFV-CPimex/main.F90   charge-preserving FV-IMEX method
|
|-- CMakeLists.txt                 method/grid selection and build rules
|-- build.sh                       command-line build helper
|-- run.slurm                      cluster-specific SLURM example
`-- README.md
```

Each method directory contains one complete main program. Shared configuration
and output code is kept under `src/common`. Files with identical module names
under `src/grids/standard` and `src/grids/half-grid` are alternative
implementations; CMake compiles exactly one grid directory into each
executable.

## Build

Run all commands from the `PICfortran` directory.

Using the helper script:

```bash
./build.sh --method TSI-Sim
./build.sh --method TSI-Sex --build-dir build-sex
./build.sh --method TSI-CPex --build-dir build-cpex
./build.sh --method TSFV-CPimex --build-dir build-tsfv
```

Using CMake directly:

```bash
cmake -S . -B build \
  -DPIC_METHOD=TSI-CPex \
  -DPIC_GRID=auto \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The executable is written to:

```text
build/bin/pic-vm.exe
```

`PIC_GRID=auto` selects the required grid automatically:

- `standard` for `TSI-Sim`, `TSI-Sex`, and `TSI-CPex`;
- `half-grid` for `TSFV-CPimex`.

CMake rejects incompatible method/grid combinations. If FFTW is installed in
a nonstandard location, pass its installation prefix explicitly:

```bash
cmake -S . -B build \
  -DPIC_METHOD=TSI-CPex \
  -DFFTW_ROOT=/path/to/fftw
```

## Run

The executable reads `config/pic.nml` by default:

```bash
mpirun -np 4 build/bin/pic-vm.exe
```

Set `PIC_CONFIG` to use another namelist:

```bash
PIC_CONFIG=config/example-two-stream.nml \
  mpirun -np 4 build/bin/pic-vm.exe

PIC_CONFIG=config/example-weibel.nml \
  mpirun -np 4 build/bin/pic-vm.exe
```

The two `example-*.nml` files are small runnable examples. They are not the
production input files used to generate manuscript figures.

## Runtime configuration

Each configuration file contains two Fortran namelists.

### `pic`

The `pic` group controls:

- the initial condition and number of particles per cell;
- the spatial, velocity, and fast-variable resolutions;
- the spatial and velocity-domain lengths;
- `epsilon`, the time step, and final time;
- the two-stream and Weibel distribution parameters;
- the wave number and initial magnetic perturbation.

The available initial-condition IDs are:

| ID | Initial condition |
|---:|---|
| 1 | counter-streaming Maxwellians |
| 2 | perturbed anisotropic Maxwellian (Weibel instability) |
| 3 | bump-on-tail distribution |

When `cfg_domain_length` is negative, the code uses one wavelength,
`2*pi/cfg_wave_number`. The number of fast-variable nodes `cfg_ntau` must be a
positive even integer.

The initial longitudinal electric field is computed from the discrete Poisson
equation. The transverse electric field initially vanishes, and the magnetic
field uses the configured sinusoidal perturbation.

### `pic_output`

The `pic_output` group controls:

- phase-space distributions `fxvx` and `fvxvy`;
- charge density, velocity moments, and field snapshots;
- mass, charge, momentum, and energy diagnostics;
- kinetic-energy components and field norms;
- physical and two-scale Poisson residuals;
- snapshot and diagnostic output intervals;
- an optional filename prefix.

Snapshot cadence is set in physical time by `cfg_snapshot_interval`.
Diagnostic cadence is set in time steps by
`cfg_diagnostic_interval_steps`.

## Output format

Output files are created below the directory from which the executable is
launched. Run different experiments in separate working directories or use a
different `cfg_output_prefix` for each case.

Scalar time series contain two columns:

```text
time  value
```

This applies to mass, charge, total energy, and Poisson residuals.

Multi-component time series contain physical components only, without a
redundant time column:

| Output | Columns |
|---|---|
| `momentum` | `P1`, `P2` |
| `kinetic_components` | `K1`, `K2` |
| `field_l2` | `||E1||_L2`, `||E2||_L2`, `||B||_L2` |

For these files, reconstruct time from the row index:

```text
time = row_index * cfg_dt * cfg_diagnostic_interval_steps
```

with `row_index=0` for the initial value.

## Method notes

### Symmetric methods

`TSI-Sim` and `TSI-Sex` use a collocated field grid and standard current
deposition. They are uniformly accurate under the well-prepared two-scale
initialization, but they do not enforce the discrete Gauss constraint exactly
during long-time integration.

### Charge-preserving exponential method

`TSI-CPex` modifies the longitudinal electric-field update using the discrete
continuity relation. The longitudinal field is initialized from the deposited
two-scale density, so the discrete Gauss constraint is satisfied initially and
propagated by the time integrator.

### Charge-preserving finite-volume method

`TSFV-CPimex` places the longitudinal electric field and current on cell faces
and the other fields on grid points. Its longitudinal current combines the
particle path current in slow time with the fast-variable current. This
construction satisfies the compatible discrete continuity equation and
preserves the staggered-grid Gauss constraint.

## Reproducibility notes

- Random-number streams depend on the MPI rank. Runs are deterministic for a
  fixed MPI decomposition, but changing the number of ranks changes the
  particle sample.
- Record the compiler, MPI implementation, FFTW version, process count, and
  complete namelist for production runs.
- The supplied SLURM script contains cluster-specific partitions, modules, and
  paths and must be adapted to another system.

## Citation and license

Add the final paper citation and DOI here when available. Before publishing the
repository, also add an explicit open-source `LICENSE` file.
