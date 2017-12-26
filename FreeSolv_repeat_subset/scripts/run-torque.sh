#!/bin/bash
#  Batch script for mpirun job on cbio cluster.
#
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=72:00:00
#
# join stdout and stderr
#PBS -j oe
#
# spool output immediately
#PBS -k oe
#
# specify queue
#PBS -q gpu
#
# nodes: number of 8-core nodes
#   ppn: how many cores per node to use (1 through 8)
#       (you are always charged for the entire node)
#PBS -l nodes=1:ppn=1:gpus=1:shared
#
# export all my environment variables to the job
##PBS -V
#
# job name (default = name of script file)
#PBS -N smirnoff

cd $PBS_O_WORKDIR
cat $PBS_GPUFILE
echo $CUDA_VISIBLE_DEVICES

source activate smirnoff

# Run the simulation. Distinguish between array
# of jobs and multiple singular jobs.
if [ -n "$PBS_ARRAYID" ]; then
  echo "Running array job $PBS_ARRAYID/${n_jobs}..."
  python run_yank.py $PBS_ARRAYID ${n_jobs}
else
  echo "Running single job ${job_id}/${n_jobs}"
  python run_yank.py ${job_id} ${n_jobs}
fi

