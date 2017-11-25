# Validation of modifications to smirnoff99Frosst

This contains code and work relating to two significant modifications of smirnoff99Frosst done at a fairly late stage in the process:
- Fixing an issue where cyclobutane was planar (inherited from parm@Frosst) which is incorrect
- Validating our new hydroxyl hydrogen parameters, which give polar hydroxyl hydrogens a small but nonzero radius and epsilon value

To validate these changes we are repeating a subset of our density, dielectric constant, and hydration free energy calculations for those compounds affected.

Here we pick out compounds using the affected parameters from the relevant test sets so that we can repeat the validation calculations on the relevant subsets.

## Procedure for determining molecules to retain
Here we loop over molecules in FreeSolv and the ThermoML density/dielectric constant validation set and see which molecules use an affected parameter from smirnoff99Frosst; those which use an affected parameter are stored to new files in a `filtered_data` directory.

## Manifest
- [`original_data`](original_data): Original data from the density/dielectric constant validation set (currently https://github.com/kyleabeauchamp/B2) and FreeSolv (github.com/mobleylab/FreeSolv) v0.52 for filtering
- [`filtered_data`](filtered_data): Filtered data from the same test sets, containing only those molecules which have modified parameters in the update.
