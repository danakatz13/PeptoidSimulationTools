#!/bin/bash

#SBATCH --job-name=amber
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

module load amber/openmpi/intel/22.00

srun sander.MPI -O -i min_test.in -p complex.parm7 -c complex.rst7 -o min1.out -r min_1.rst7
