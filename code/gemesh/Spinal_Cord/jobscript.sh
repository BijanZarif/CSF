#!/bin/bash
# Job name:
#SBATCH --job-name=test_run
#
# Project:
#SBATCH --account=nn9279k
# Wall clock limit:
#SBATCH --time='120:00:00'
#
# Max memory usage per task:
#SBATCH --mem-per-cpu=1000M
#
# Number of tasks (cores):
#SBATCH --nodes=1 --ntasks=16
#SBATCH --hint=compute_bound
#SBATCH --cpus-per-task=4

##SBATCH --partition=long
#SBATCH --output=test_run.out

## Set up job environment
source /cluster/bin/jobsetup


module load gcc/4.9.2
module load openmpi.gnu/1.8.4
#source ~oyvinev/intro/hashstack/fenics-1.5.0.abel.gnu.conf
source ~oyvinev/fenics1.5/fenics1.5

# Expand pythonpath with locally installed packages
export PYTHONPATH=$PYTHONPATH:$HOME/.local/lib/python2.7/site-packages/

# Define what to do when job is finished (or crashes)
cleanup "mkdir -p /work/users/vegarvi/intro"
cleanup "cp -r $SCRATCH/RESULTS /work/users/vegarvi/intro/results"

# Copy necessary files to $SCRATCH
cp script.py $SCRATCH/
cp FSI_Solver.py $SCRATCH/
cp longer_top_spinal.xml $SCRATCH/

# Enter $SCRATCH and run job
cd $SCRATCH
mpirun --bind-to none python script.py
