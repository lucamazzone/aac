#!/bin/bash -l
#SBATCH --job-name=attempting
#SBATCH --time=09:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --output=aout_a.out
#SBATCH --error=aerr_a.err
#SBATCH --constraint=mc

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1
module load daint-gpu

srun python attempt.py

