#!/bin/bash

#SBATCH --job-name=production
#SBATCH --gres=gpu:1
#SBATCH --mem=20GB
#SBATCH --time=24:00:00

module load amber/openmpi/intel/22.00

pmemd.cuda -O -i md.in -p complex.parm7 -c heat.rst7 -o md1.mdout -x md1.nc -r md1.rst7
