#!/bin/bash -l

#SBATCH --ntasks=40
#SBATCH -N 2
#SBATCH --time=10:15:00
#SBATCH --output=mapinvert_output.txt
#SBATCH --error=mapinvert_error.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=luca.g.mazzone@gmail.com

####################

mpiexec -np $SLURM_NTASKS python using.py