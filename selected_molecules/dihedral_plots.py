#!/bin/env python

import mdtraj
import numpy as np
import os
import matplotlib
from pylab import *

# Which molecules to analyze
molnames = ['testcase_smiles1', 'testcase_name3']
# Where to write output graphs
outdir = 'dihedral_plots'
if not os.path.isdir(outdir): os.mkdir(outdir)


# Store which dihedrals to analyze by molecule name
dihedrals = {}
# testcase_smiles1:
# - proper involving nitrogen ring - oxygen-carbon-carbon-carbon(adjacent to nitrogen); measures whether nitrogen-containing ring is flipping
dihedrals['testcase_smiles1'] = [[2, 13, 11, 8]]
# - proper involving sulfur -- carbon-carbon-(bridgehead)-carbon-sulfur; measures whether sulfur-containing ring is flipping.
dihedrals['testcase_smiles1'].append( [6, 12, 14, 3] )

# testcase_name3
# - torsion involving four aromatic carbons within central ring (all bridgeheads); buckles in GAFF/GAFF2
dihedrals['testcase_name3'] = [[26, 28, 29, 27 ]]



# Loop over molecules and compute dihedrals
for molname in molnames:
    # Loop over trajectory types
    for ttype in ['', '_gaff', '_gaff2']:
        # Specify trajectory name to load
        trajfile = 'trajectories/%s%s.nc' % (molname,ttype)
        # Load a PDB file specifying the topology and the trajectory
        traj = mdtraj.load(trajfile, top='minimized/%s_smirff.pdb' % molname)

        for idx, atomset in enumerate(dihedrals[molname]):
            # Compute dihedrals
            dvals = mdtraj.compute_dihedrals(traj, np.array(atomset, ndmin=2))

            figure()
            hist(dvals, 36)
            xlabel('dihedral angle (radians')
            ylabel('Counts')
            xticks([-pi/2, -pi/4, 0, pi/4, pi/2], [r'$-\frac{\pi}{2}$', r'$-\frac{\pi}{4}$', '$0$', r'$-\frac{\pi}{4}$', r'$\frac{\pi}{2}$'])
            #ylabel
            savefig( os.path.join(outdir, '%s%s_%s.pdf' % (molname, ttype, idx)))

            # TODO: Make axis ticks/labels invlolve pi
            # use fixed axis range?
