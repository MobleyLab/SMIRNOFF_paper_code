#!/usr/local/bin/env python

import os
import sys
import glob
import json

from yank.yamlbuild import YamlBuilder


MOLECULES_DONE_FILEPATH = os.path.join('..', 'results', 'molecules_done.json')


def print_and_flush(msg):
    """On the cluster, printing may not always appear if not flushed."""
    print(msg)
    sys.stdout.flush()


def read_status():
    if os.path.exists(MOLECULES_DONE_FILEPATH):
        print_and_flush('Node {}: Resuming from {}'.format(job_id, MOLECULES_DONE_FILEPATH))
        with open(MOLECULES_DONE_FILEPATH, 'r') as f:
            molecules_done = set(json.load(f))
    else:
        molecules_done = set()
    return molecules_done


if __name__ == '__main__':
    job_id = int(sys.argv[1])
    n_job = int(sys.argv[2])
    assert 0 < job_id <= n_job

    # Load cached status calculations.
    molecules_done = read_status()

    # Get all experiments to run and sort them so that order is deterministic.
    yank_scripts_filepaths = glob.glob(os.path.join('..', 'yank', '*.yaml'))
    yank_scripts_filepaths = sorted(yank_scripts_filepaths)

    # Get all YANK YAML scripts.
    for i, yank_script_filepath in enumerate(yank_scripts_filepaths):
        # Check if the job is assigned to this script
        # and/or if we have already completed this.
        if (i  % n_job != job_id - 1 or yank_script_filepath in molecules_done):
            print_and_flush('Node {}: Skipping {}'.format(job_id, yank_script_filepath))
            continue

        # Run YANK.
        print_and_flush('Node {}: Running {}'.format(job_id, yank_script_filepath))
        yaml_builder = YamlBuilder(yank_script_filepath)
        yaml_builder.run_experiments()

        # Update completed molecules. We need to load this every
        # time because parallel nodes may have updated this. This
        # is subject to race conditions, but in the rare case it
        # fails, the script should fix it automatically in the next
        # run.
        molecules_done = read_status()
        molecules_done.add(yank_script_filepath)
        with open(MOLECULES_DONE_FILEPATH, 'w') as f:
            json.dump(list(molecules_done), f)
