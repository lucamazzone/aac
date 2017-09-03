#!/bin/bash -l
#SBATCH --job-name=using
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=2
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --output=out_asg.out
#SBATCH --error=err_asg.err
#SBATCH --constraint=mc


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun python using.py


