#!/bin/bash -l
#SBATCH --job-name=attempting
#SBATCH --time=07:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --account=s555
#SBATCH --output=out_asg.out
#SBATCH --error=err_asg.err


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1
module load daint-gpu

srun python using.py

