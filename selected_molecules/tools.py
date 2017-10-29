#!/bin/env/python

from simtk import unit
import numpy as np
import simtk.openmm as mm
from simtk.openmm import app
from simtk.openmm import Platform
from openeye.oechem import *

#==============================================================================
# UTILITY FUNCTIONS
#==============================================================================

def setPositionsInOEMol(molecule, positions):
    """Set the positions in an OEMol using a position array with units from simtk.unit, i.e. from OpenMM. Atoms must have same order.

    Arguments:
    ---------
    molecule : OEMol
        OpenEye molecule
    positions : Nx3 array
        Unit-bearing via simtk.unit Nx3 array of coordinates
    """
    if molecule.NumAtoms() != len(positions): raise ValueError("Number of atoms in molecule does not match length of position array.")
    pos_unitless = positions/unit.angstroms

    coordlist = []
    for idx in range(len(pos_unitless)):
        for j in range(3):
            coordlist.append( pos_unitless[idx][j])
    molecule.SetCoords(OEFloatArray(coordlist))



def extractPositionsFromOEMol(molecule):
    """Get the positions from an OEMol and return in a position array with units via simtk.unit, i.e. foramtted for OpenMM.
    Adapted from choderalab/openmoltools test function extractPositionsFromOEMOL

    Arguments:
    ----------
    molecule : OEMol
        OpenEye molecule

    Returns:
    --------
    positions : Nx3 array
        Unit-bearing via simtk.unit Nx3 array of coordinates
    """

    positions = unit.Quantity(np.zeros([molecule.NumAtoms(), 3], np.float32), unit.angstroms)
    coords = molecule.GetCoords()
    for index in range(molecule.NumAtoms()):
        positions[index,:] = unit.Quantity(coords[index], unit.angstroms)
    return positions

def minimizeOpenMM(Topology, System, Positions):
    """
    Minimize molecule with specified topology, system, and positions
       using OpenMM. Return the positions of the optimized moelcule.

    Parameters
    ----------
    Topology
    System
    Positions

    Returns
    -------
    minimized positions

    """

    # need to create integrator but don't think it's used
    integrator = mm.LangevinIntegrator(
            300.0 * unit.kelvin,
            1.0 / unit.picosecond,
            2.0 * unit.femtosecond)

    simulation = app.Simulation(Topology, System, integrator)
    simulation.context.setPositions(Positions)
    simulation.minimizeEnergy()

    Topology.positions = simulation.context.getState(getPositions=True).getPositions()

    return Topology.positions

def reformat_oemol_coordinates( oemol ):
    """
    Take an oemol with multiple conformers and return a conformer_positions numpy array which is an Nx3xM array
    where N is the number of atoms, each has three coordinates, and M is the number of conformers.

    Parameters
    ----------
    oemol : OpenEye oemol, multiconformer
        Multi-conformer molecule to work with

    Returns
    ----------
    conformer_positions : numpy array
        Nx3xM array of simtk.unit.Quantity of dimension (natoms,3, M) where M is the number of quantity, with units of length
        """

    #Get positions for use below
    numconfs = oemol.NumConfs()
    coordinates = oemol.GetCoords()
    natoms=len(coordinates)

    conformer_positions = np.zeros([natoms,3, numconfs], np.float32)
    for (idx,conf) in enumerate(oemol.GetConfs()):
        for index in range(natoms):
            coordinates=conf.GetCoords()

            (x,y,z) = coordinates[index]
            conformer_positions[index,0, idx] = x
            conformer_positions[index,1, idx] = y
            conformer_positions[index,2, idx] = z
    conformer_positions = unit.Quantity(conformer_positions, unit.angstroms)
    return conformer_positions
