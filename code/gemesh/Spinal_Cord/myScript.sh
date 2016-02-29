#!/bin/bash

# Job name:
#SBATCH --job-name=test_dt_E
#
# Project:
#SBATCH --account=nn9279k
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=4G

## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

## Copy input files to the work directory:
cp tests_dt_E/FSI_solver.py $SCRATCH
cp tests_dt_E/script.py $SCRATCH
cp tests_dt_E/longer_top_spinal.xml $SCRATCH

## Make sure the results are copied back to the submit directory (see Work Directory below):
chkfile RESULTS

## Do some work:
cd $SCRATCH
YourCommands
