# SMIRNOFF paper code
Code/tools relating to the initial SMIRNOFF format paper.

## Manifest
* [`FreeSolv/`][FreeSolv]: Tools/scripts relating to FreeSolv hydration free energy calculations as reported in the paper; see additional README file there
* [`FreeSolv_repeat_subset/`][FreeSolv_repeat_subset]: Tools/scripts for repeating a subset of the FreeSolv hydration free energy calculations with an updated forcefield (changes to cyclobutane and hydroxyl radii); see additional README file there. Mirrors `FreeSolv` but for a subset of molecules selected in `modification_validation`
* [`selected_molecules`](selected_molecules): Scripts/analysis relating to gas phase simulations (and torsional analysis) of selected molecules as reported in the paper.
* [`jupyter_notebooks`](jupyter_notebooks): Jupyter notebooks utilized in the paper.
* [`hydrogen_radii`](hydrogen_radii): Work relating to adding nonzero LJ radii/well-depths to hydroxyl hydrogens to fix issues with AMBER family force fields.
* [`modification_validation`](modification_validation): Code/tools used to validate modification of smirnoff99Frosst on a subset of density/dielectric/hydration calculations after fixes to cyclobutane parameters and adding hydroxyl radii.
