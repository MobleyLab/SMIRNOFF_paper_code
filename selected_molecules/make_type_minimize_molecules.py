#!/bin/env python

from openeye.oechem import *
from openeye.oeomega import *
from openeye.oeiupac import *
import os
import numpy as np
from openmoltools.openeye import *

import simtk.openmm as mm
from simtk.openmm import app
from simtk.openmm import Platform
from simtk import unit

from tools import * #Helper stuff in this directory
#==============================================================================
# MOLECULES TO MAKE FROM NAMES
#==============================================================================

#Biphenyl bridge problem cases
iupaclist = ['biphenyl', '1,2-diphenylbenzene', '1,2,3-triphenylbenzene', '1,2,3,4-tetraphenylbenzene']

#==============================================================================
# MOLECULS TO MAKE FROM SMILES
#==============================================================================
# Rotatable bond problem case/challenging atom typing: methylenepenta-1,4-diene
smileslist = [ 'C(=C)C(=C)C(=C)' ]

# Bridgehead examples from Bayly, May 2017
smileslist += ['c1coc(c3c[nH]cc3)c1c2sccc2']

#==============================================================================
# GENERATE MOLECULES
#==============================================================================

mols = {}
ct=0

# For molecules in iupac list, create molecule, ID, and store to dict
for name in iupaclist:
    mol = OEMol()
    OEParseIUPACName(mol, name)
    smi = OECreateIsoSmiString(mol)
    mol.SetTitle(name)
    molid = 'testcase_name%s' % ct
    mols[molid] = [ mol, name, smi ]
    ct+=1

# For molecules in SMILES list, create molecule, name, and store to dict
ct=0
for smi in smileslist:
    mol = OEMol()
    OEParseSmiles(mol, smi)
    mol.SetTitle(name)
    name = OECreateIUPACName(mol)
    OETriposAtomNames(mol)
    molid = 'testcase_smiles%s' % ct
    mols[molid] = [ mol, name, smi ]
    ct+=1


#==============================================================================
# GENERATE CONFORMERS AND CHARGES
#==============================================================================
# Conformers here will be generated with OpenEye's Omega

# Generate OEMol
mol = OEMol()

# Initialize omega
omega = OEOmega()
omega.SetStrictStereo(False) #Choose random stereoisomer for those with unspecified stereochemistry
omega.SetMaxConfs(1)

# Loop over input molecules, generate conformers and write
if not os.path.isdir('outmol2_tripos'): os.mkdir('outmol2_tripos')
if not os.path.isdir('outmol2_gaff'): os.mkdir('outmol2_gaff')
if not os.path.isdir('outmol2_gaff2'): os.mkdir('outmol2_gaff2')
for molid in mols:
    #Retrieve molecule, generate conf
    mol = mols[molid][0]
    omega(mol)
    #Assign charges
    mol = get_charges(mol)
    #Re-store molecule
    mols[molid][0] = mol
    ofs = oemolostream('outmol2_tripos/%s.mol2' % molid)
    OEWriteMolecule(ofs, mol)
    ofs.close()

    # Also make gaff mol2 file
    os.system('antechamber -i %s/%s.mol2 -fi mol2 -o %s/%s.mol2 -fo mol2 -at gaff' % ('outmol2_tripos', molid, 'outmol2_gaff', molid))
    # Make frcmod file for gaff
    os.system('parmchk2 -i %s/%s.mol2 -f mol2 -o %s/%s.frcmod' % ('outmol2_gaff', molid, 'outmol2_gaff', molid))
    # Also make gaff2 mol2 file
    os.system('antechamber -i %s/%s.mol2 -fi mol2 -o %s/%s.mol2 -fo mol2 -at gaff2' % ('outmol2_tripos', molid, 'outmol2_gaff2', molid))
    # Also make gaff2 frcmod
    os.system('parmchk2 -i %s/%s.mol2 -f mol2 -s gaff2 -o %s/%s.frcmod' % ('outmol2_gaff2', molid, 'outmol2_gaff2', molid))


#==============================================================================
# GENERATE AMBER PARAMETER FILES
#==============================================================================
# This utilizes the Chodera lab's `openmoltools` as a driver for tleap to parameterize the GAFF and GAFF2 cases and generate AMBER format input files.


from openmoltools.amber import *

if not os.path.isdir('prmtop_gaff'): os.mkdir('prmtop_gaff')
if not os.path.isdir('prmtop_gaff2'): os.mkdir('prmtop_gaff2')

for molid in mols:
    print("Setting up %s for AMBER" % molid)
    # Gaff
    run_tleap(molid, 'outmol2_gaff/%s.mol2' % molid, 'outmol2_gaff/%s.frcmod' % molid, prmtop_filename = 'prmtop_gaff/%s.prmtop' % molid, inpcrd_filename = 'prmtop_gaff/%s.crd' % molid)
    # Gaff 2
    run_tleap(molid, 'outmol2_gaff2/%s.mol2' % molid, 'outmol2_gaff2/%s.frcmod' % molid, prmtop_filename = 'prmtop_gaff2/%s.prmtop' % molid, inpcrd_filename = 'prmtop_gaff2/%s.crd' % molid, leaprc='leaprc.gaff2')


#==============================================================================
# ENERGY MINIMIZE DIFFERENT CASES
#==============================================================================
# This uses OpenMM to energy minimize the different cases, reading the AMBER parameter files generated above in the GAFF/GAFF2 cases, and assigning SMIRNOFF parameters from smirnoff99frosst.ffxml in the SMIRNOFF case.


# Prep/file bookkeeping - GAFF case
gaff_indirs = ['prmtop_gaff', 'prmtop_gaff2']
outdir = 'minimized'
if not os.path.isdir(outdir): os.mkdir(outdir)

# Prep force field for smirnoff case
import parmed
from openforcefield.typing.engines.smirnoff import *
from parmed.openmm import topsystem
ff = ForceField('forcefield/smirnoff99Frosst.ffxml')

topologies = {}
systems = {}

# Actually do minimizations
print("Minimizing systems:")
for molid in mols:
    print("   Minimizing %s with SMIRFF..." % molid)
    # Prep for minimization using SMIRFF
    mol = mols[molid][0]
    topology = generateTopologyFromOEMol(mol)
    topologies[molid] = topology
    cpositions = reformat_oemol_coordinates(mol)
    system = ff.createSystem(topology, [mol], verbose=False)
    systems[molid] = system

    # Put into parmed
    parm = topsystem.load_topology( topology, system=None, xyz=None, box=None)
    # Minimize
    parm.positions =minimizeOpenMM(topology, system, cpositions[:,:,0])
    # Write valid mol2 file, bypasssing ParmEd which won't have bond order info
    ofs = oemolostream( os.path.join(outdir, '%s_smirff.mol2' % molid))
    setPositionsInOEMol(mol, parm.positions)
    OEWriteMolecule(ofs, mol)
    ofs.close()

    # Write pdb file
    pdbf = parmed.formats.PDBFile
    pdbf.write(parm, os.path.join(outdir, '%s_smirff.pdb' % molid) )

    # DO MINIMIZATIONS USING GAFF PARAMETERS (GAFF AND GAFF2)
    # Start w GAFF
    print("      with GAFF...")
    parm = parmed.load_file('prmtop_gaff/%s.prmtop' % molid, 'prmtop_gaff/%s.crd' % molid)
    Topology = parm.topology
    Positions = parm.positions
    System = parm.createSystem(nonbondedMethod=app.NoCutoff)
    parm.positions = minimizeOpenMM( Topology, System, Positions)
    # Use PDB file; reasoning explained above
    pdbf = parmed.formats.PDBFile
    pdbf.write(parm, os.path.join(outdir, '%s_gaff.pdb' % molid) )


    # GAFF2
    print("      with GAFF2...")
    parm = parmed.load_file('prmtop_gaff2/%s.prmtop' % molid, 'prmtop_gaff2/%s.crd' % molid)
    Topology = parm.topology
    Positions = parm.positions
    System = parm.createSystem(nonbondedMethod=app.NoCutoff)
    parm.positions = minimizeOpenMM( Topology, System, Positions)
    # Write PDB file
    pdbf = parmed.formats.PDBFile
    pdbf.write(parm, os.path.join(outdir, '%s_gaff2.pdb' % molid) )
