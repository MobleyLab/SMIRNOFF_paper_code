#!/usr/local/bin/env python

import os
import glob
import json

import mdtraj as md
from simtk import openmm, unit
from openeye import oechem
from openmoltools import packmol
from openforcefield.utils import get_data_filename
from openforcefield.typing.engines import smirnoff


def create_pdb_files(molecule_smiles, molecule_mol2_filepath, water_mol2_filepath,
                     vacuum_pdb_filepath, solvated_pdb_filepath):
    """Create vacuum and solvated molecule box and save in PDB format.

    The function uses packmol to solvate the molecule in 1309 waters.

    Parameters
    ----------
    molecule_smiles : str
        SMILES string for the molecule.
    molecule_mol2_filepath : str
        Path to the mol2 file describing the molecule to solvate.
    water_mol2_filepath : str
        Path to the mol2 file describing the water molecule geometry.
    vacuum_pdb_filepath : str
        Output path for the PDB file describing the molecule in vacuum.
    solvated_pdb_filepath : str
        Output path for the PDB file of the solvated molecule.

    """
    n_monomers = [1, 1309]  # molecule, waters
    mol_traj = md.load(molecule_mol2_filepath)
    water_traj = md.load(water_mol2_filepath)

    # Export molecule in PDB format for vacuum calculation.
    mol_traj.save_pdb(vacuum_pdb_filepath)

    # Estimate box size based on water density.
    box_size = packmol.approximate_volume_by_density([molecule_smiles, 'O'], n_monomers)

    # Solvate molecule and export PDB.
    packed_traj = packmol.pack_box([mol_traj, water_traj], n_molecules_list=[1, 1309],
                                   box_size=box_size)
    packed_traj.save_pdb(solvated_pdb_filepath)


def load_oe_graph_mol(mol2_filepath):
    """Load a single OpenEye molecule from a mol2 file."""
    # OpenEye flavors to use when loading mol2 files.
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield

    ifs = oechem.oemolistream(mol2_filepath)
    ifs.SetFlavor(oechem.OEFormat_MOL2, flavor)
    mol = oechem.OEGraphMol()

    molecules = []
    while oechem.OEReadMolecule(ifs, mol):
        oechem.OETriposAtomNames(mol)
        molecules.append(oechem.OEGraphMol(mol))

    # The script assumes we have only 1 molecule per mol2 file.
    assert len(molecules) == 1
    return molecules[0]


def create_xml_system_files(molecule_mol2_filepath, tip3p_mol2_filepath,
                            vacuum_pdb_filepath, solvated_pdb_filepath,
                            vacuum_xml_filepath, solvated_xml_filepath):
    """Create vacuum and solvated OpenMM systems and save in XML format.

    Parameters
    ----------
    molecule_mol2_filepath : str
        Path to the mol2 file describing the molecule.
    tip3p_mol2_filepath : str
        Path to the mol2 file describing the water molecule with TIP3P charges.
    vacuum_pdb_filepath : str
        Path to the PDB file for the molecule in vacuum.
    solvated_pdb_filepath : str
        Path to the PDB file for the solvated molecule.
    vacuum_xml_filepath : str
        Output path of the vacuum system.
    solvated_xml_filepath : str
        Output path of the solvated system.

    """
    # Load TIP3P water into OpenEye molecule.
    water_oe_mol = load_oe_graph_mol(tip3p_mol2_filepath)

    # Load solvated molecule into OpeneEye molecule.
    molecule_oe_mol = load_oe_graph_mol(molecule_mol2_filepath)

    # Load forcefield.
    # TODO add HBonds constraint through createSystem when openforcefield#32 is implemented.
    ff = smirnoff.ForceField('hydrogen_radii.ffxml', 'hbonds.ffxml', 'forcefield/tip3p.ffxml')

    # Load vacuum and solvated PDBFiles.
    pdb_file_vacuum = openmm.app.PDBFile(vacuum_pdb_filepath)
    pdb_file_solvated = openmm.app.PDBFile(solvated_pdb_filepath)

    # Create vacuum and solvated system.
    system_vacuum = ff.createSystem(pdb_file_vacuum.topology, molecules=[molecule_oe_mol],
                                    nonbondedMethod=smirnoff.NoCutoff)
    # TODO add HBonds constraints here when openforcefield#32 is implemented.
    system_solvated = ff.createSystem(pdb_file_solvated.topology, molecules=[molecule_oe_mol, water_oe_mol],
                                      nonbondedMethod=smirnoff.PME, nonbondedCutoff=1.1*unit.nanometer,
                                      ewaldErrorTolerance=1e-4)  #, constraints=smirnoff.HBonds)

    # Fix switching function.
    # TODO remove this when openforcefield#31 is fixed
    for force in system_solvated.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.setUseSwitchingFunction(True)
            force.setSwitchingDistance(1.0*unit.nanometer)

    # Save XML files.
    for system, xml_filepath in zip([system_vacuum, system_solvated],
                                    [vacuum_xml_filepath, solvated_xml_filepath]):
        system_serialized = openmm.XmlSerializer.serialize(system)
        with open(xml_filepath, 'w') as f:
            f.write(system_serialized)


def create_yank_input_files():
    """Creates XML, PDB and YAML script files to run YANK calculations on a subset of the FreeSolv set.

    Starting from mol2 files, the script uses packmol to solvate the molecules
    in TIP3P waters and produces PDB files of the vacuum and solvated boxes
    together with OpenMM XML System files parametrized with smirnoff99Frosst.
    A YANK script is generated for each system following the unified FreeSolv
    protocol.

    """
    database_filepath = os.path.join('..', 'database.json')
    smirnoff_openmm_dir = os.path.join('..', 'openmmfiles')
    pdb_dir = os.path.join('..', 'pdbfiles')
    water_mol2_filepath = get_data_filename(os.path.join('systems', 'monomers', 'tip3p_water.mol2'))

    # Create directories that don't exist.
    for directory_path in [smirnoff_openmm_dir, pdb_dir]:
        if not os.path.exists(directory_path):
            os.makedirs(directory_path)
    # Read in database.
    with open(database_filepath, 'r') as f:
        database = json.load(f)

    # Find all molecules.
    mol2_files_pattern = os.path.join('..', 'mol2files_sybyl', '*.mol2')
    tripos_filenames = list(glob.glob(mol2_files_pattern))

    # Remove everything with names which are NOT in the list of those we need to repeat.
    file = open( '../../modification_validation/filtered_data/FreeSolv/database_filtered.txt', 'r')
    database_text = file.readlines()
    file.close()
    retain_mols = []
    for line in database_text:
        #If it's not a comment, parse it and add to list
        if line[0] != '#':
            retain_mols.append(line.split(';')[0])

    fnm = []
    for filename in tripos_filenames:
        molid = os.path.basename(filename).replace('.mol2','')
        if molid in retain_mols:
            fnm.append(filename)
    tripos_filenames = fnm

    # Create YANK input files.
    for i, tripos_filename in enumerate(tripos_filenames):
        molecule_id, _ = os.path.splitext(os.path.basename(tripos_filename))  # Remove path and extension.
        print('\rSetting up {} ({}/{})     '.format(molecule_id, i+1, len(tripos_filenames)))

        # Output file paths.
        vacuum_filename = molecule_id + '_vacuum'
        solvated_filename = molecule_id + '_solvated'
        vacuum_pdb_filepath = os.path.join(pdb_dir, vacuum_filename + '.pdb')
        solvated_pdb_filepath = os.path.join(pdb_dir, solvated_filename + '.pdb')
        vacuum_smirnoff_openmm_filepath = os.path.join(smirnoff_openmm_dir, vacuum_filename + '.xml')
        solvated_smirnoff_openmm_filepath = os.path.join(smirnoff_openmm_dir, solvated_filename + '.xml')

        # Export solvated and vacuum molecule in PDB format.
        create_pdb_files(database[molecule_id]['smiles'], tripos_filename, water_mol2_filepath,
                         vacuum_pdb_filepath, solvated_pdb_filepath)

        # Create OpenMM XML file parametrized with SMIRNOFF.
        create_xml_system_files(tripos_filename, water_mol2_filepath,
                                vacuum_pdb_filepath, solvated_pdb_filepath,
                                vacuum_smirnoff_openmm_filepath, solvated_smirnoff_openmm_filepath)

if __name__ == '__main__':
    create_yank_input_files()
