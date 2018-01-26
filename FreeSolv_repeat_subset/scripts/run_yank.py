#!/usr/local/bin/env python

import os
import sys
import glob
import json

from yank.yamlbuild import YamlBuilder


RESULTS_DIR = os.path.join('..', 'results')
MOLECULES_DONE_FILEPATH = os.path.join(RESULTS_DIR, 'molecules_done.json')


def print_and_flush(msg):
    """On the cluster, printing may not always appear if not flushed."""
    print(msg)
    sys.stdout.flush()


def read_status():
    if os.path.exists(MOLECULES_DONE_FILEPATH):
        print_and_flush('Node {}: Resuming from {}'.format(job_id, MOLECULES_DONE_FILEPATH))
        with open(MOLECULES_DONE_FILEPATH, 'r') as f:
            molecules_data = set(json.load(f))
    else:
        molecules_data = set()
    return molecules_data


def update_status(new_element):
    """Update status with the new piece of data.

    We need to load this every time because parallel nodes may have
    updated this. This is subject to race conditions, but in the rare
    case it fails, the script should fix it automatically in the next
    run.
    """
    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)

    molecules_data = read_status()
    molecules_data.add(new_element)
    with open(MOLECULES_DONE_FILEPATH, 'w') as f:
        json.dump(list(molecules_data), f)


def run_yank(job_id, n_jobs):
    openmm_system_dir = os.path.join('..', 'openmmfiles')
    pdb_dir = os.path.join('..', 'pdbfiles')
    yank_script_template_filepath = 'yank_template.yaml'

    # Read in YANK template script.
    with open(yank_script_template_filepath, 'r') as f:
        script_template = f.read()

    # Load cached status calculations.
    molecules_done = read_status()

    # Find all molecules to run.
    molecules_files_pattern = os.path.join(pdb_dir, '*_vacuum.pdb')
    molecule_ids = [os.path.basename(molecule_file)[:-11]
                    for molecule_file in glob.glob(molecules_files_pattern)]

    # Sort molecules so that parallel nodes won't make the same calculation.
    molecule_ids = sorted(molecule_ids)

    # Create YANK input files.
    for i, molecule_id in enumerate(molecule_ids):

        # Check if the job is assigned to this script and/or if we
        # have already completed this.
        if (i % n_jobs != job_id - 1 or
                    molecule_id in molecules_done):
            print_and_flush('Node {}: Skipping {}'.format(job_id, molecule_id))
            continue

        # Output file paths.
        vacuum_filename = molecule_id + '_vacuum'
        solvated_filename = molecule_id + '_solvated'
        vacuum_pdb_filepath = os.path.join(pdb_dir, vacuum_filename + '.pdb')
        solvated_pdb_filepath = os.path.join(pdb_dir, solvated_filename + '.pdb')
        vacuum_xml_filepath = os.path.join(openmm_system_dir, vacuum_filename + '.xml')
        solvated_xml_filepath = os.path.join(openmm_system_dir, solvated_filename + '.xml')

        # Create yank script.
        phase1_path = str([solvated_xml_filepath, solvated_pdb_filepath])
        phase2_path = str([vacuum_xml_filepath, vacuum_pdb_filepath])
        script = script_template.format(experiment_dir=molecule_id,
                                        phase1_path=phase1_path, phase2_path=phase2_path)

        # Run YANK.
        print_and_flush('Node {}: Running {}'.format(job_id, molecule_id))
        yaml_builder = YamlBuilder(script)
        yaml_builder.run_experiments()

        # Update completed molecules.
        update_status(molecule_id)


if __name__ == '__main__':
    job_id = int(sys.argv[1])
    n_jobs = int(sys.argv[2])
    assert 0 < job_id <= n_jobs
    run_yank(job_id=job_id, n_jobs=n_jobs)
