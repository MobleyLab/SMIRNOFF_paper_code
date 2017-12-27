#!/bin/bash
N_JOBS=10
# for i in 30 31 32 33
for ((i=1; i<=$N_JOBS; i++ ))
do
  qsub -v job_id=$i,n_jobs=$N_JOBS run-torque.sh
done
