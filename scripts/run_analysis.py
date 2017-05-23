#!/usr/local/bin/env python

import os
import json
import logging

import numpy as np


# A validation test fails when its Z-score exceeds this threshold.
MAX_Z_SCORE = 6
ANALYSIS_FILEPATH = os.path.join('..', 'results', 'analysis_done.json')

# Set logger verbosity level.
logging.basicConfig(level=logging.DEBUG)


def analyze_directory(experiment_dir):
    """Return free energy and error of a single experiment.

    Parameters
    ----------
    experiment_dir : str
        The path to the directory storing the nc files.

    Return
    ------
    DeltaF : simtk.unit.Quantity
        Difference in free energy between the end states in kcal/mol.
    dDeltaF: simtk.unit.Quantity
        Statistical error of the free energy estimate.

    """
    import yaml
    from simtk import unit
    from yank.analyze import get_analyzer

    analysis_script_filepath = os.path.join(experiment_dir, 'analysis.yaml')

    # Load sign of alchemical phases.
    with open(analysis_script_filepath, 'r') as f:
        analysis_script = yaml.load(f)

     # Generate analysis object.
    analysis = {}
    for phase_name, sign in analysis_script:
        phase_path = os.path.join(experiment_dir, phase_name + '.nc')
        analyzer = get_analyzer(phase_path)
        analysis[phase_name] = analyzer.analyze_phase()
        kT = analyzer.kT

    # Compute free energy.
    DeltaF = 0.0
    dDeltaF = 0.0
    for phase_name, sign in analysis_script:
        DeltaF -= sign * (analysis[phase_name]['DeltaF'] + analysis[phase_name]['DeltaF_standard_state_correction'])
        dDeltaF += analysis[phase_name]['dDeltaF']**2
    dDeltaF = np.sqrt(dDeltaF)

    # Convert from kT units to kcal/mol
    unit_conversion = kT / unit.kilocalories_per_mole
    return DeltaF * unit_conversion, dDeltaF * unit_conversion


# @mpi.on_single_node(rank=0)
def print_analysis(experiment_name, expected_free_energy, obtained_free_energy):
    """Print the results of the analysis.

    Parameters
    ----------
    experiment_name : str
        The name of the experiment to print.
    expected_free_energy : tuple of float
        The expected pair (DeltaF, dDeltaF) in kT.
    obtained_free_energy : tuple of float
        The pair (DeltaF, dDeltaF) in kT from the calculation.

    """
    expected_DeltaF, expected_dDeltaF = expected_free_energy
    obtained_DeltaF, obtained_dDeltaF = obtained_free_energy

    # Determine if test has passed.
    z_score = (obtained_DeltaF - expected_DeltaF) / expected_dDeltaF
    test_passed = abs(z_score) < MAX_Z_SCORE

    # Print results.
    print('{}: {}\n'
          '\texpected: {:.3f}, {:.3f} kcal/mol\n'
          '\tobtained: {:.3f}, {:.3f} kcal/mol\n'
          '\tZ-score: {:.3f}'.format(experiment_name, 'OK' if test_passed else 'FAIL',
                                     expected_DeltaF, expected_dDeltaF,
                                     obtained_DeltaF, obtained_dDeltaF,
                                     z_score))


# @mpi.on_single_node(rank=0)
def save_analysis(analysis, filepath):
    with open(filepath, 'w') as f:
        json.dump(analysis, f)


def read_freesolv(molecule_ids):
    # Import FreeSolv database free energies.
    with open(os.path.join('..', 'database.json'), 'r') as f:
        freesolv_database = json.load(f)

    # Isolate calculated free energies and uncertainties.
    freesolv_database = {molecule_id: [freesolv_database[molecule_id]['calc'],
                                       freesolv_database[molecule_id]['d_calc'],
                                       freesolv_database[molecule_id]['expt'],
                                       freesolv_database[molecule_id]['d_expt']]
                         for molecule_id in molecule_ids}
    return freesolv_database


def run_analysis():
    """Run analysis on all validation tests."""

    # Load completed analysis.
    if os.path.exists(ANALYSIS_FILEPATH):
        with open(ANALYSIS_FILEPATH, 'r') as f:
            analysis_done = json.load(f)
    else:
        analysis_done = dict()

    # Load completed calculations.
    with open('molecules_done.json', 'r') as f:
        molecules_done = json.load(f)

    # Molecules done is now a list of paths to the YAML files. Convert to molecules ids.
    for i in range(len(molecules_done)):
        molecule_id = os.path.basename(molecules_done[i])
        molecule_id = os.path.splitext(molecule_id)[0]
        molecules_done[i] = molecule_id

    # Isolate calculated free energies and uncertainties.
    freesolv_database = read_freesolv(set(molecules_done))

    # Filter molecules with those that we have already analyzed.
    molecules_to_analyze = {molecule_id for molecule_id in molecules_done
                            if molecule_id not in analysis_done}

    # # Distribute analysis of remaining experiments across nodes.
    # from yank import mpi
    # experiment_dirs = [os.path.join('..', 'experiments', molecule_id)
    #                    for molecule_id in molecules_to_analyze]
    # free_energies = mpi.distribute(analyze_directory, experiment_dirs, send_results_to='all')
    #
    # # Store computed free energies.
    # for i, molecule_id in enumerate(molecules_to_analyze):
    #     analysis_done[molecule_id] = free_energies[i]
    # save_analysis(analysis_done, ANALYSIS_FILEPATH)

    # Run analysis and store results.
    for molecule_id in molecules_to_analyze:
        experiment_dir = os.path.join('..', 'experiments', molecule_id)
        free_energy = analyze_directory(experiment_dir)
        analysis_done[molecule_id] = free_energy
        save_analysis(analysis_done, ANALYSIS_FILEPATH)

    # Print comparisons.
    for molecule_id, free_energy in analysis_done.items():
        print_analysis(molecule_id, freesolv_database[molecule_id][:2], free_energy)


def plot_correlation():
    from matplotlib import pyplot as plt

    # Read analysis and freesolv database.
    with open(ANALYSIS_FILEPATH, 'r') as f:
        analysis = json.load(f)
    molecule_ids = list(analysis.keys())
    freesolv_database = read_freesolv(analysis)

    # Prepare data for plotting.
    gaff_f = [freesolv_database[molecule_id][0] for molecule_id in molecule_ids]
    gaff_df = [freesolv_database[molecule_id][1] for molecule_id in molecule_ids]
    expt_f = [freesolv_database[molecule_id][2] for molecule_id in molecule_ids]
    expt_df = [freesolv_database[molecule_id][3] for molecule_id in molecule_ids]
    smirnoff_f = [analysis[molecule_id][0] for molecule_id in molecule_ids]
    smirnoff_df = [analysis[molecule_id][1] for molecule_id in molecule_ids]

    # Plot scatter.
    fig, (ax1, ax2) = plt.subplots(ncols=2, sharex=True, figsize=(13, 6))

    # Add horizontal line and labels.
    for ax, name, ref_f, ref_df in zip([ax1, ax2], ['GAFF', 'Exp'], [gaff_f, expt_f], [gaff_df, expt_df]):
        ax.errorbar(ref_f, smirnoff_f, xerr=ref_df, yerr=smirnoff_df, fmt='o')
        ax.set(adjustable='box-forced', aspect='equal')
        x = np.array(ax.get_xlim())
        y = np.array(ax.get_ylim())
        ax.fill_between(x, y - 2, y + 2, alpha=0.4)
        ax.fill_between(x, y - 1, y + 1, alpha=0.4)
        ax.plot(x, y, ls='--', c='black')
        ax.set_xlabel('{} $\Delta F$ [kcal/mol]'.format(name))
        ax.set_ylabel('SMIRNOFF $\Delta F$ [kcal/mol]')

    plt.show()
    # plt.savefig('results.pdf')


if __name__ == '__main__':
    run_analysis()
    # plot_correlation()
