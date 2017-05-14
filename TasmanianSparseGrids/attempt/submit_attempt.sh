#!/bin/bash -l

#SBATCH --ntasks=32
#SBATCH --time=20:15:00
#SBATCH --output=mapinvert_output.txt
#SBATCH --error=mapinvert_error.err
#SBATCH --partition=hiprioq
#SBATCH --nodelist=node12



#SBATCH --mail-type=ALL
#SBATCH --mail-user=luca.g.mazzone@gmail.com

####################

mpiexec -np $SLURM_NTASKS python using.py