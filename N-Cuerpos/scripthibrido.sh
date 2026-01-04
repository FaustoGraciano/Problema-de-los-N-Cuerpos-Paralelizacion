#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH -o directorioSalida/output.txt
#SBATCH -e directorioSalida/errors.txt
mpirun --bind-to none hibrido $1 $2 $3 $4