#!/bin/bash

#SBATCH --job-name=heating
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --time=72:00:00

module load amber/openmpi/intel/22.00

srun sander.MPI -O -i heat.in -p complex.parm7 -c min1.rst -o heat.mdout -x heat.nc -r heat.rst7
