#!/bin/bash

#SBATCH --job-name=amber
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --time=72:00:00

module load amber/openmpi/intel/22.00

mpirun -np 48 sander.MPI -ng 12 -groupfile remd.groupfile

`
