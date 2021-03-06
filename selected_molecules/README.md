# Comparison of sample molecules across GAFF, GAFF2, and smirnoff99Frosst

This provides some basic scripts used for the SMIRNOFF paper for comparison of gas phase energy minimizations/simulations for some sample challenging molecules in GAFF, GAFF2, and smirnoff99Frosst.
These are not production-level tools, but were used in running some simple tests as reported in the paper.

## Manifest

### Scripts
- `make_type_minimize_molecules.py`: Prepares molecules specified by IUPAC names or SMILES strings, generating conformers, assigning parameters, writing intermediate files, and energy minimizing in the gas phase with GAFF, GAFF2, and smirnoff99Frosst. Output files will be named based on the origin -- naming schemes are `testcase_[identifier]NN*` where the identifier is either `smiles` for molecules generated from SMILES strings or `name` for molecules generated from IUPAC name, and NN is the number of molecule (counting from zero), as provided in the script. Outputs are described further below. 
- `run_molecule.py`: Very simple MD simulation script which takes a selected molecule as output by `make_type_minimize_molecules.py` and runs gas phase MD simulations of the selected molecule, storing trajectories to NetCDF files (`.nc`). Trajectory length (currently, and as in plots) 5e6 steps of 2 fs each, so 10 ns.
- `dihedral_plots.py`: Plots selected dihedral angles for the two problem cases shown in th the paper, the 1,2,3,4-tetraphenyhlbenzene case and the bridgehead problem case.
- `tools.py`: Helper tools used by the above two scripts, such as for working with positions from OpenEye molecules/exchanging these with positions from OpenMM topologies.

### Outputs
- `outmol2_gaff`, `outmol2_gaff2`, `outmol2_tripos`: Output files of molecules after initial creation and conformer generation with OpenEye's Omega; `outmol2_tripos` provides files directly out of the OpenEye toolkits, and `outmol2_gaff` and `outmol2_gaff2` provide files after processing via antechamber to assign GAFF/GAFF2 atom types. The latter two directories also include AMBER frcmod files provided by `parmchk2`
- `prmtop_gaff` and `prmtop_gaff2`: AMBER format input files for the GAFF/GAFF2 cases as generated by `tleap` as driven by the Chodera lab's `openmoltools` Python module.
- `minimized`: Energy minimized molecules (mol2 format for SMIRNOFF case, PDB format for GAFF/GAFF2 case) resulting from the above.
- `trajectories`: Output NetCDF trajectories as generated by `run_molecule.py`
- `pdbfdiles`: Output PDF format version of trajectories as output by `nc_to_pdb.py`
- `dihedral_plots`: Plots of selected dihedrals (as further described in `dihedral_plots.py`) for the problem cases noted in the paper.

### Molecule names
In this particular case, molecules are named by origin (created from name or SMILES). The output molecules shown in the paper have the following partial file names:
- methylenepenta-1,4-diene: testcase_smiles0
- biphenyl: testcase_name0
- 1,2-diphenylbenzene: testcase_name1 
- 1,2,3,4-tetraphenylbenzene: testcase_name3
- bridgehead problem case (IUPAC name not very useful): test_smiles1
The "ring problem case" shown in the figure was taken from the GAFF paper and not constructed here because it was already identified in the GAFF paper as a failure case.

## Prerequisites
Use of these tools will require an installation of the following:
- AmberTools or `ambermini` (from Omnia)
- OpenEye's Python toolkits and a valid OpenEye license
- OpenMM
- `openforcefield`
- ParmEd
- `openmoltools` from the Chodera lab
- `mdtraj`
- plus associated prerequisites (NumPy, etc.)
Many of these can be conda installed from the Omnia conda channel, and `conda install -c omnia openforcefield` plus `conda install -c omnia openmoltools` will probably provide all or almost all of the required functionality, with the exception of needing your own installation of the OpenEye toolkits.

## Usage

If you would like to use the tools provided here to reproduce the work done in the paper, do the following:
- Run `make_type_minimize_molecules.py` to make, type, and minimize the molecules selected in the script (which cover the cases highlighted in the paper)
- Run `run_molecule.py` on (NOTE WHICH CASES) via `python run_molecule.py ID` to run and store trajectories
- Run `dihedral_plots.py` to do dihedral angle plots for `testcase_name3` (bridgehead problem case from the paper) and `testcase_smiles1` (1,2,3,4-tetraphenhylbenzene) as reported in the supporting information
- Run `nc_to_pdb.py` to generate PDB format versions of the NetCDF files, suitable for visualization in PyMol or elsewhere; these were used in creating visualizations of representative structures for the tetraphenyl and bridgehead cases shown in the paper.

If you would like to use this to examine your own test molecules, you would want to edit `make_type_minimize_molecules.py` to add generation of your molecules (by adding molecules to the IUPAC or SMILES lists), and note the IDs which will result for the molecules you choose to generate. Then run `run_molecule.py` on these IDs to generate trajectories, run `nc_to_pdb.py` to convert trajectories to PDB format if desired, etc. If you wish to examine specific dihedral angles, identify the indices of atoms involved in the torsions you are interested in and modify `dihedral_plots.py` to plot those dihedral angles.
