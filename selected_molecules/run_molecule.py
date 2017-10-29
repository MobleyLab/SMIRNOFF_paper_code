#!/bin/env python

import time
import simtk.openmm as mm
from simtk.openmm import app
from simtk.openmm import Platform
from simtk.unit import *
import numpy as np
from mdtraj.reporters import NetCDFReporter
from openforcefield.typing.engines.smirnoff import *
import os
from tools import *
import sys

#Define what molecule to work on, and a few simulation parameters
#molname = 'testcase_smiles1'
molname = sys.argv[1]
mol_filename =  molname+'.mol2'
time_step = 2 #Femtoseconds
temperature = 300 #kelvin
friction = 1 # per picosecond
num_steps = 5000000
trj_freq = 1000 #steps
data_freq = 1000 #steps

# Output trajectory directory
outdir = 'trajectories'
if not os.path.isdir(outdir): os.mkdir(outdir)

# Load OEMol
ifs = oechem.oemolistream( os.path.join('minimized/%s_smirff.mol2' % molname))
mol = oechem.OEGraphMol()
#ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
oechem.OEReadMolecule(ifs, mol )
oechem.OETriposAtomNames(mol)

# Get positions
coordinates = mol.GetCoords()
natoms = len(coordinates)
positions = np.zeros([natoms,3], np.float32)
for index in range(natoms):
    (x,y,z) = coordinates[index]
    positions[index,0] = x
    positions[index,1] = y
    positions[index,2] = z
positions = Quantity(positions, unit.angstroms)

# Load forcefield
forcefield = ForceField('forcefield/smirnoff99Frosst.ffxml')

#Define system
topology = generateTopologyFromOEMol(mol)
system = forcefield.createSystem(topology, [mol])

#Do simulation - SMIRFF
integrator = mm.LangevinIntegrator(temperature*kelvin, friction/picoseconds, time_step*femtoseconds)
#platform = mm.Platform.getPlatformByName('Reference')
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature*kelvin)
netcdf_reporter = NetCDFReporter(os.path.join(outdir, molname+'.nc'), trj_freq)
simulation.reporters.append(netcdf_reporter)
simulation.reporters.append(app.StateDataReporter('data.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True))

print("Starting simulation for %s" % molname)
start = time.clock()
simulation.step(num_steps)
end = time.clock()

print("Elapsed time %.2f seconds" % (end-start))
netcdf_reporter.close()
print("Done!")

#Do simulation - GAFF
import parmed
parm = parmed.load_file('prmtop_gaff/%s.prmtop' % molname, 'prmtop_gaff/%s.crd' % molname)
topology, positions = parm.topology, parm.positions
system = parm.createSystem(nonbondedMethod=app.NoCutoff)
integrator = mm.LangevinIntegrator(temperature*kelvin, friction/picoseconds, time_step*femtoseconds)
#platform = mm.Platform.getPlatformByName('Reference')
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature*kelvin)
netcdf_reporter = NetCDFReporter(os.path.join(outdir, molname+'_gaff.nc'), trj_freq)
simulation.reporters.append(netcdf_reporter)
simulation.reporters.append(app.StateDataReporter('data_gaff.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True))

print("Starting GAFF simulation")
start = time.clock()
simulation.step(num_steps)
end = time.clock()

print("Elapsed time %.2f seconds" % (end-start))
netcdf_reporter.close()
print("Done!")

#Do simulation - GAFF2
import parmed
parm = parmed.load_file('prmtop_gaff2/%s.prmtop' % molname, 'prmtop_gaff2/%s.crd' % molname)
topology, positions = parm.topology, parm.positions
system = parm.createSystem(nonbondedMethod=app.NoCutoff)
integrator = mm.LangevinIntegrator(temperature*kelvin, friction/picoseconds, time_step*femtoseconds)
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature*kelvin)
netcdf_reporter = NetCDFReporter(os.path.join(outdir, molname+'_gaff2.nc'), trj_freq)
simulation.reporters.append(netcdf_reporter)
simulation.reporters.append(app.StateDataReporter('data_gaff2.csv', data_freq, step=True, potentialEnergy=True, temperature=True, density=True))

print("Starting GAFF2 simulation")
start = time.clock()
simulation.step(num_steps)
end = time.clock()

print("Elapsed time %.2f seconds" % (end-start))
