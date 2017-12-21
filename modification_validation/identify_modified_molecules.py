#!/bin/env python
from openeye.oechem import *
from openforcefield.typing.engines.smirnoff.forcefield import *
import os
import pandas as pd

# List parameters which were modified here; molecules containing these
# will be retrieved
modified_smirks = [ '[#1:1]-[#8]', '[#6r4:1]-;@[#6r4:2]-;@[#6r4:3]', '!#1:1]-[#6r4:2]-;!@[!#1:3]', '[!#1:1]-[#6r4:2]-;!@[#1:3]']

# output directory
outdir = 'filtered_data'
if not os.path.isdir(outdir): os.mkdir(outdir)
# input directory
indir = 'original_data'

# Set up forcefield
ff = ForceField('forcefield/smirnoff99Frosst.ffxml')

# Define utility functionality to tell us if a particular molecule uses a particular SMIRKS
def molecule_uses_parameter_bysmirks(mol, smirks):
    """Check whether a provided OEMol uses a parameter specified by its SMIRKS string.
    Arguments:
    ----------
        mol: OEMol of molecule
        smirks: SMIRKS pattern in question; all force types will be checked

    Returns:
    --------
        Boolean value as to whether that SMIRKS pattern is used for a parameter applied to that molecule.
    """
    # Label molecules
    labels = ff.labelMolecules( [ mol] )
    # Retrieve just this molecule's labels
    labels = labels[0]

    # Extract SMIRKS patterns used
    used = False
    for ftype in labels:
        for fterm in labels[ftype]:
            if fterm[2]==smirks:
                #print("FOUND")
                #print(fterm)
                used = True
                return True
    return used


# Process FreeSolv case
file = open( os.path.join(indir, 'FreeSolv/database.txt'), 'r')
database_text = file.readlines()
file.close()
targetdir = os.path.join(outdir, 'FreeSolv')
if not os.path.isdir(targetdir): os.mkdir(targetdir)
outfile = open( os.path.join(targetdir, 'database_filtered.txt'), 'w')

for line in database_text:
    #If it's not a comment, parse it and skip ahead if it's not a molecule we want
    if line[0] != '#':
        # Extract SMILES
        smi = line.split(';')[1]
        # Build OEMol for this SMILES
        mol = OEMol()
        OESmilesToMol(mol, smi)
        OEAddExplicitHydrogens(mol)
        # Check if this SMILES is affected by parameter change
        used = False
        for smirks in modified_smirks:
            if molecule_uses_parameter_bysmirks(mol, smirks):
                used = True
        # Don't store if not affected, otherwise store (below)
        if not used:
            continue

    # Write to output if it's a comment or a molecule we want to keep
    outfile.write(line)

outfile.close()


# Process density/dielectrics
# Easiest way to load data from these files is with pandas.read_csv but I'm not
# that proficient with pandas yet. There may be a better way to do this, but
# here's what I've got.

input_databases = [os.path.join(indir, 'density_dielectric', 'data_dielectric.csv'),
                    os.path.join(indir, 'density_dielectric', 'full_filtered_data.csv')]
source_dfs = [ pd.read_csv(d) for d in input_databases]

# Loop over dataframes. For each, build a set of SMILES, check to see which
# are affected by our new parameters, and dump a new dataframe
# of affected instances
for (idx, df) in enumerate(source_dfs):
    # Make set of all SMILES
    all_smiles = { v for k, v in df['smiles'].items()}
    affected_smiles = []

    # Build OEMol for each SMILES and check
    for smi in all_smiles:
        mol = OEMol()
        OESmilesToMol(mol, smi)
        OEAddExplicitHydrogens(mol)
        # Check if this SMILES is affected by parameter change
        used = False
        for smirks in modified_smirks:
            if molecule_uses_parameter_bysmirks(mol, smirks):
                affected_smiles.append(smi)
                break

    # Dump dataframe of affected cases
    # Create dataframe
    new_df = df.loc[df['smiles'].isin(affected_smiles)]
    # Prep output file
    targetdir=os.path.join(outdir, 'density_dielectric')
    if not os.path.isdir(targetdir): os.mkdir(targetdir)
    filename = os.path.basename( input_databases[idx].replace('.csv', '_filtered.csv'))
    # Write output
    new_df.to_csv( os.path.join( targetdir, filename))
