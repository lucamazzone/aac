#!/bin/bash -l
#SBATCH --job-name=main
#SBATCH --mail-user=luca.g.mazzone@gmail.com
#SBATCH --time=04:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --output=out_asg.out
#SBATCH --error=err_asg.err


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1
module load daint-gpu


srun ./mainfile.exec

