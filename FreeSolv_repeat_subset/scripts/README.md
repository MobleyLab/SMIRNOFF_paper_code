# YANK and analysis scripts
The scripts generally require the [`database.json`](https://github.com/MobleyLab/FreeSolv/blob/master/database.json)
file in the FreeSolv repo to be stored in the upper folder.

## Manifest
* `create_input_files.py`: create PDB and OpenMM XML files to run the simulation.
* `yank_template.yaml`: basic YAML script specifying the protocol of the calculations for all molecules.
* `run_yank.py`: run YANK calculations.
* `run-torque.sh`: bash script to execute `run_yank.py` on a Torque cluster.
* `run_analysis.py`: script to run the analysis of the YANK calculations; in this particular case it analyzes just the subset of calculations re-run in this directory (i.e. for hydroxyls and cyclobutyl compounds); also extended to analyze the SAME subset of compounds but from the original directory for cross-comparison.
* `run_analysis_fullset.py`: modified version of run_analysis.py which pulls the remaining compounds of the set from the `FreeSolv` directory and thus analyzes the full set (original calculations for the non-hydroxyl/cyclobutyl compounds, new calculations for the hydroxyl/cyclobutyl compounds).
* `hbonds.ffxml`: SMIRNOFF file to constrain hydrogen bonds.
