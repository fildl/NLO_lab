# NEMO Tides Experiment

> [!IMPORTANT]
> **University Project Disclaimer**: This repository contains code for a university project. It is **not** intended for public replication or general use. The simulations are configured to run on a specific remote cluster environment and **will not work locally** without significant modification.
>
> **Usable Content**: The only part of this repository that is intended to be portable and potentially usable by others is the **data analysis** section (`postprocessing/`).

This repository contains the configuration and analysis scripts for a NEMO (Nucleus for European Modelling of the Ocean) experiment simulating tidal/Kelvin wave propagation in an idealized basin.

## Project Structure

- **`TIDES/`**: Contains the NEMO source code modifications and experiment configuration.
    - `MY_SRC/`: User-defined Fortran source files (e.g., `usrdef_*.F90`) defining the domain, initial state, and boundary conditions.
    - `EXP00/`: Run directory containing namelists (`namelist_cfg`) and execution scripts.
- **`postprocessing/`**: Contains Python scripts for analyzing and visualizing the simulation results.
    - `experiments/`: Directory for experiment-specific data/configs.
    - `plot_kelvin_wave.py`: Script to generate Hovmöller diagrams and SSH snapshots.

## Experiment Description

The experiment is designed to simulate a Kelvin wave propagating in an **idealized basin representing the Adriatic Sea**. Key configurations include:
- **Domain**: Idealized rectangular basin with flat bottom (defined in `usrdef_zgr.F90`), mimicking the Adriatic geometry.
- **Forcing**: Initial SSH perturbation (Gaussian bump) or tidal forcing (defined in `usrdef_istate.F90` or `usrdef_sbc.F90`).
- **Physics**: Linear free surface, no wind stress (unless modified).

## Getting Started

### Prerequisites

- **NEMO** (v4.0 or compatible version installed).
- **Fortran Compiler** (e.g., gfortran, ifort).
- **Python 3** with standard scientific stack (`numpy`, `matplotlib`, `netCDF4` or `xarray`).

### Running the Simulation

> [!NOTE]
> The simulation configuration is tailored for a specific remote cluster. The following steps are for reference and documentation purposes only and are not expected to work in a local environment.

1.  Navigate to the experiment directory:
    ```bash
    cd TIDES2/EXP00
    ```
2.  (Optional) Clean and rebuild if source code in `MY_SRC` was modified:
    ```bash
    ../makenemo -m <ARCH> -n TIDES2 -r NEMO_REF -d "MY_SRC"
    ```
    *Note: Adjust `<ARCH>` to your machine's architecture file.*
3.  Run the simulation:
    ```bash
    ./run_nemo_v2.sh
    ```
    *Or run the executable directly if no script is preferred.*

### Analysis

The analysis scripts found in `postprocessing/` are independent of the cluster environment and can be run locally to visualize the results (provided you have the output data).

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

[University Project]
This material is part of a university course project and is provided for educational demonstration purposes only. It is not intended for commercial use or replication in production environments.
