# NEMO Tides Experiment

This repository contains the configuration and analysis scripts for a NEMO (Nucleus for European Modelling of the Ocean) experiment simulating tidal/Kelvin wave propagation in an idealized basin.

## Project Structure

- **`TIDES/`**: Contains the NEMO source code modifications and experiment configuration.
    - `MY_SRC/`: User-defined Fortran source files (e.g., `usrdef_*.F90`) defining the domain, initial state, and boundary conditions.
    - `EXP00/`: Run directory containing namelists (`namelist_cfg`) and execution scripts.
- **`postprocessing/`**: Contains Python scripts for analyzing and visualizing the simulation results.
    - `experiments/`: Directory for experiment-specific data/configs.
    - `plot_kelvin_wave.py`: Script to generate Hovmöller diagrams and SSH snapshots.

## Experiment Description

The experiment is designed to simulate a Kelvin wave propagating in a rectangular basin. Key configurations include:
- **Domain**: Idealized rectangular basin with flat bottom (defined in `usrdef_zgr.F90`).
- **Forcing**: Initial SSH perturbation (Gaussian bump) or tidal forcing (defined in `usrdef_istate.F90` or `usrdef_sbc.F90`).
- **Physics**: Linear free surface, no wind stress (unless modified).

## Getting Started

### Prerequisites

- **NEMO** (v4.0 or compatible version installed).
- **Fortran Compiler** (e.g., gfortran, ifort).
- **Python 3** with standard scientific stack (`numpy`, `matplotlib`, `netCDF4` or `xarray`).

### Running the Simulation

1.  Navigate to the experiment directory:
    ```bash
    cd TIDES/EXP00
    ```
2.  (Optional) Clean and rebuild if source code in `MY_SRC` was modified:
    ```bash
    ../makenemo -m <ARCH> -n TIDES -r NEMO_REF -d "MY_SRC"
    ```
    *Note: Adjust `<ARCH>` to your machine's architecture file.*
3.  Run the simulation:
    ```bash
    ./run_nemo_v2.sh
    ```
    *Or run the executable directly if no script is preferred.*

### Analysis

1.  Navigate to the postprocessing directory:
    ```bash
    cd postprocessing
    ```
2.  Run the plotting script:
    ```bash
    python plot_kelvin_wave.py
    ```
3.  Outputs (plots) will be generated in the same directory (Note: `.png` files are git-ignored by default).

## License

[Specify License if applicable, e.g., MIT, Student Project]
