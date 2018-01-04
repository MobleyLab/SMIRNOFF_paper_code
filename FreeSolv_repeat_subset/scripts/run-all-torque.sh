#!/bin/bash
#PBS -n manager
#PBS -q home-gibbs
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00

N_JOBS=50
# for i in 30 31 32 33
for ((i=1; i<=$N_JOBS; i++ ))
do
  qsub -v job_id=$i,n_jobs=$N_JOBS run-torque.sh
done
