# Density/dielectric constant benchmarks

## Overview

As part of this study, we ran a benchmark on densities/dielectric constants of a set of molecules from NIST's ThermoML, as previously reported in this study from the Chodera lab.
Code for that study is available at [LiquidBenchmark](https://github.com/choderalab/LiquidBenchmark) on the Chodera Lab GitHub site; however, dependencies have gone out of date and the code needed cleanup, so author Kyle Beauchamp helped get the code working again and made it available at [https://github.com/kyleabeauchamp/SolutionFFBench](https://github.com/kyleabeauchamp/SolutionFFBench).
We also worked with him to extend it to be compatible with the SMIRNOFF format.
Here, we used versions 0.0.1 and 0.0.2 of this updated benchmark code to run and analyze our density and dielectric constant calculations for GAFF (v0.0.1) and SMIRNOFF (v0.0.2).
Those releases of Kyle Beauchamp's code are archived here.
Additionally, we provide run settings used when conducting the simulations.


## Manifest
- [`B2-0.0.1.tar.gz`](B2-0.0.1.tar.gz): Version 0.0.1 of benchmark code as used for GAFF benchmarks. Repo subsequently renamed to `SolutionFFBench`
- [`B2-0.0.2.tar.gz`](B2-0.0.2.tar.gz): Version 0.0.2 of benchmark code as used for SMIRNOFF benchmarks. Repo subsequently renamed to `SolutionFFBench`
- [`results_GAFF`](results_GAFF): Results of GAFF calculations (tables and plots)
- [`results_newSMIRNOFF`](results_newSMIRNOFF): Results of SMIRNOFF calculations (tables and plots)


## Settings
Default settings deposited in the `B2-*.tar.gz` files are not those used for running production level simulations here.
Here, we used (in `density_simulation_parameters.py`):
- `N_STEPS = 1000000`
- `N_EQUIL_STEPS = 10000000`
- `OUTPUT_FREQUENCY_EQUIL = 100000`
- `OUTPUT_FREQUENCY = 5000`
- `OUTPUT_DATA_FREQUENCY = 250`
- `STD_ERROR_TOLERANCE = 0.0002` (g/mL)
The forcefield name was varied from GAFF to SMIRNOFF for output files and to deterermine the input forcefield.
