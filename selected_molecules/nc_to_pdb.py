#!/bin env python

# Convert NetCDF files for cases output here to PDB format for easier visualization across different viewers

import mdtraj
import numpy as np
import os
import glob

# Create output directory if not present
if not os.path.isdir('pdbfiles'): os.mkdir('pdbfiles')

# Get list of trajectory file names to convert
trajfiles_gaff2 = glob.glob('trajectories/*_gaff2.nc')
traj_prefixes = []
for name in trajfiles_gaff2:
    newname= name.replace('_gaff2.nc', '')
    traj_prefixes.append(newname)

# Loop over trajectories and convert
for ttype in ['', '_gaff', '_gaff2']:
    for tprefix in traj_prefixes:
        molname = os.path.basename(tprefix)
        # Load relevant NetCDF trajectory
        trajfile = 'trajectories/%s%s.nc' % (molname,ttype)
        # We require a topology file which has the chemical contents of the system, so use a PDB file stored previously
        traj = mdtraj.load(trajfile, top='minimized/%s_smirff.pdb' % molname)

        # Align trajectory to frame 0 and write to PDB for viewing
        traj.superpose(traj[0])
        traj.save(os.path.join('pdbfiles','%s%s.pdb' % (molname,ttype)))


