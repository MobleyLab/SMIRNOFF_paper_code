#!/usr/local/bin/env python

import os
import json
import logging

import numpy as np
from scipy.stats import linregress


# A validation test fails when its Z-score exceeds this threshold.
MAX_Z_SCORE = 6
ANALYSIS_FILEPATH = os.path.join('..', 'results', 'analysis_done.json')
FULL_ANALYSIS_FILEPATH = os.path.join('../../', 'FreeSolv/results', 'analysis_done.json')
MOLECULES_DONE_FILEPATH = os.path.join('..', 'results', 'molecules_done.json')
FULL_MOLECULES_DONE_FILEPATH = os.path.join('../..', 'FreeSolv/results', 'molecules_done.json')

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
        n_samples = len(analyzer.reporter._storage_dict['analysis'].dimensions['iteration']) - 1
        assert n_samples == 5000, '{}'.format(n_samples)

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
                                       freesolv_database[molecule_id]['d_expt'],
                                       freesolv_database[molecule_id]['iupac']]
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


    # Load completed old analysis and add any compounds we DON'T have to dict
    if os.path.exists(FULL_ANALYSIS_FILEPATH):
        with open(FULL_ANALYSIS_FILEPATH, 'r') as f:
            old_analysis_done = json.load(f)
            # Copy in entries
            for key in old_analysis_done:
                if not key in analysis_done:
                    analysis_done[key] = old_analysis_done[key]

    # Load completed calculations.
    with open(MOLECULES_DONE_FILEPATH, 'r') as f:
        molecules_done = json.load(f)
    print(len(molecules_done))
    with open(FULL_MOLECULES_DONE_FILEPATH, 'r') as f:
        old_molecules_done = json.load(f)
        for cid in old_molecules_done:
            if not cid in molecules_done:
                molecules_done.append(cid)

    # Remove a duplicate compound which was removed in FreeSolv v0.52
    molecules_done.remove('mobley_4689084')

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
    interrupted_run = {}
    for molecule_id in molecules_to_analyze:
        experiment_dir = os.path.join('..', 'experiments', molecule_id)
        try:
            free_energy = analyze_directory(experiment_dir)
        except AssertionError as e:
            interrupted_run[molecule_id] = str(e)
            print('\nFound interrupted calculation: molecule {}\n'.format(molecule_id))
        #except FileNotFoundError:
            # LOAD OLD DATA
        else:
            analysis_done[molecule_id] = free_energy
            save_analysis(analysis_done, ANALYSIS_FILEPATH)
    # Remove molecule we don't want to analyze
    analysis_done.pop('mobley_4689084')

    # Print comparisons.
    for molecule_id, free_energy in analysis_done.items():
        print_analysis(molecule_id, freesolv_database[molecule_id][:2], free_energy)

    # Print interrupted molecules.
    print('The calculations of the following molecules have been interrupted.')
    for molecule_id, n_samples in interrupted_run.items():
        print(molecule_id, n_samples)

    # Print list version of interrupted molecules for convenience.
    print()
    print('["' + '", "'.join(interrupted_run.keys()) + '"]')


def compute_sample_statistics(x, y):
    slope, intercept, r_value, p_value, stderr = linregress(x, y)
    diff = np.array(x) - np.array(y)
    avg_err = diff.mean()
    rms_err = np.sqrt((diff**2).mean())
    return np.array([r_value, avg_err, rms_err])


def compute_bootstrap_statistics(x, y, percentile=0.95, n_bootstrap_samples=5000):
    """Compute R, mean error and RMSE for the data with bootstrap confidence interval."""
    x = np.array(x)
    y = np.array(y)
    statistics = compute_sample_statistics(x, y)

    # Generate bootstrap statistics variations.
    statistics_samples_diff = np.zeros((n_bootstrap_samples, len(statistics)))
    for i in range(n_bootstrap_samples):
        bootstrap_sample_indices = np.random.randint(low=0, high=len(x), size=len(x))
        x_bootstrap = x[bootstrap_sample_indices]
        y_bootstrap = y[bootstrap_sample_indices]
        statistics_samples_diff[i] = compute_sample_statistics(x_bootstrap, y_bootstrap)

    # Compute confidence intervals.
    stat_bound_id = int(np.floor(n_bootstrap_samples * (1 - percentile) / 2)) - 1
    statistics_confidence_intervals = []
    for i, stat_samples_diff in enumerate(statistics_samples_diff.T):
        stat_samples_diff.sort()
        stat_lower_bound = stat_samples_diff[stat_bound_id]
        stat_higher_bound = stat_samples_diff[-stat_bound_id+1]
        statistics_confidence_intervals.append([statistics[i], (stat_lower_bound, stat_higher_bound)])

    return statistics_confidence_intervals


def plot_correlation(print_outliers=False):
    from matplotlib import cm, pyplot as plt

    plt.rcParams['font.size'] = 6

    # Read analysis and freesolv database.
    with open(ANALYSIS_FILEPATH, 'r') as f:
        analysis = json.load(f)

    with open(FULL_ANALYSIS_FILEPATH, 'r') as f:
        old_analysis = json.load(f)
        # Copy in entries
        for key in old_analysis:
            if not key in analysis:
                analysis[key] = old_analysis[key]
    analysis.pop('mobley_4689084') #Compound removed from FreeSolv v0.52

    molecule_ids = list(analysis.keys())
    freesolv_database = read_freesolv(analysis)

    # Prepare data for plotting.
    gaff_f = [freesolv_database[molecule_id][0] for molecule_id in molecule_ids]
    gaff_df = [freesolv_database[molecule_id][1] for molecule_id in molecule_ids]
    expt_f = [freesolv_database[molecule_id][2] for molecule_id in molecule_ids]
    expt_df = [freesolv_database[molecule_id][3] for molecule_id in molecule_ids]
    smirnoff_f = [analysis[molecule_id][0] for molecule_id in molecule_ids]
    smirnoff_df = [analysis[molecule_id][1] for molecule_id in molecule_ids]

    # Print outliers.
    if print_outliers:
        error_threshold = 1.0  # kcal/mol
        statistical_error_theshold = 0.2  # kcal/mol
        print('# molecules: ', len(molecule_ids))
        print('# id; abs(DDf) [kcal/mol]; statistical_error [kcal/mol]; iupac name')
        for i, molecule_id in enumerate(molecule_ids):
            diff = abs(smirnoff_f[i] - gaff_f[i])
            if diff > error_threshold or smirnoff_df[i] > statistical_error_theshold:
                print('{}; {:.3f} kcal/mol; {:.3f} kcal/mol; {}'.format(molecule_id, diff, smirnoff_df[i],
                                                                        freesolv_database[molecule_id][4]))

    # Obtain colors from color map.
    colormap = cm.get_cmap('viridis')
    colors = [colormap(x) for x in [0.25, 0.5, 0.85]]

    # Plot scatter.
    fig, subplot_axes = plt.subplots(ncols=3, sharex=True, sharey=True, figsize=(7.5, 2.8))
    ax1, ax2, ax3 = subplot_axes

    gaff_vs_smirnoff = [ax1, ('GAFF', 'SMIRNOFF'), (gaff_f, smirnoff_f), (gaff_df, smirnoff_df)]
    expt_vs_smirnoff = [ax2, ('Exp', 'SMIRNOFF'), (expt_f, smirnoff_f), (expt_df, smirnoff_df)]
    expt_vs_gaff = [ax3, ('Exp', 'GAFF'), (expt_f, gaff_f), (expt_df, gaff_df)]

    # Correlation plot.
    x_limits = np.array([np.inf, -np.inf])
    for i, (ax, labels, f, df) in enumerate([gaff_vs_smirnoff, expt_vs_smirnoff, expt_vs_gaff]):
        ax.errorbar(f[0], f[1], xerr=df[0], yerr=df[1], fmt='o', markersize='1.5', color='black',  #color=colors[1],
                    alpha=0.85, elinewidth=0.75, ecolor='black', capsize=2, capthick=0.5)
        ax.set(adjustable='box-forced', aspect='equal')
        ax.set_xlabel('{} $\Delta F$ [kcal/mol]'.format(labels[0]))
        ax.set_ylabel('{} $\Delta F$ [kcal/mol]'.format(labels[1]))

        # Compute statistics and bootstrap confidence interval.
        bootstrap_statistics = compute_bootstrap_statistics(f[0], f[1])
        r_value, r_value_interval = bootstrap_statistics[0]
        avg_err, avg_err_interval = bootstrap_statistics[1]
        rms_err, rms_err_interval = bootstrap_statistics[2]
        if i == 0:
            statistics_msg = ("R$^2$:          {:.3f} [ {:.3f},  {:.3f}]\n"
                              "MD:   {:.3f} [{:.3f}, {:.3f}] kcal/mol\n"
                              "RMS diff: {:.3f} [ {:.3f},  {:.3f}] kcal/mol")
        else:
            statistics_msg = ("R$^2$:      {:.3f} [ {:.3f},  {:.3f}]\n"
                              "ME:  {:.3f} [{:.3f}, {:.3f}] kcal/mol\n"
                              "RMSE: {:.3f} [ {:.3f},  {:.3f}] kcal/mol")
        statistics_msg = statistics_msg.format(
            r_value**2, r_value_interval[0]**2, r_value_interval[1]**2,
            avg_err, avg_err_interval[0], avg_err_interval[1],
            rms_err, rms_err_interval[0], rms_err_interval[1])
        ax.text(-27.0, 2.2, statistics_msg)

        # Remove spines.
        for spine_name in ['top', 'bottom', 'right', 'left']:
            ax.spines[spine_name].set_visible(False)
        # Remove ticks.
        ax.tick_params(which='both', top=False, bottom=False, left=False, right=False,
                       labeltop=False, labelbottom=True, labelleft=True, labelright=False)
        # Add grid.
        ax.grid(axis='both', linestyle="--", lw=0.5, color="black", alpha=0.2)

        ax_limits = ax.get_xlim()
        x_limits[0] = x_limits[0] if x_limits[0] < ax_limits[0] else ax_limits[0]
        x_limits[1] = x_limits[1] if x_limits[1] > ax_limits[1] else ax_limits[1]

    # Add diagonal line now that we know the total size of the plots.
    for ax, _, _, _ in [gaff_vs_smirnoff, expt_vs_smirnoff, expt_vs_gaff]:
        ax.fill_between(x_limits, x_limits + 1, x_limits + 2, alpha=0.4, color=colors[0])
        ax.fill_between(x_limits, x_limits - 2, x_limits - 1, alpha=0.4, color=colors[0])
        ax.fill_between(x_limits, x_limits - 1, x_limits + 1, alpha=0.4, color=colors[2])
        ax.plot(x_limits, x_limits, ls='--', c='black', alpha=0.8, lw=0.7)
        ax.set_xlim(*x_limits)
        ax.set_ylim(*x_limits)


    fig.tight_layout()
    # plt.show()
    plt.savefig('../results/correlation_plots_fullset.pdf')


if __name__ == '__main__':
    run_analysis()
    plot_correlation()
    #plot_convergence()
