#!/bin/bash

#SBATCH --job-name=production
#SBATCH --gres=gpu:1
#SBATCH --mem=20GB
#SBATCH --time=24:00:00

module load amber/openmpi/intel/22.00

pmemd.cuda -O -i md-eq-nvt-5ns.in -o equil-nvt.out -p design_2.parm7 -c md_eq-npt.rst7 -r md_eq-nvt.rst7 -x equil-nvt.nc -ref md_eq-npt.rst7

