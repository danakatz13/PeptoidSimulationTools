#!/bin/bash

#SBATCH --job-name=production
#SBATCH --gres=gpu:1
#SBATCH --mem=20GB
#SBATCH --time=24:00:00

module load amber/openmpi/intel/22.00

pmemd.cuda -O -i md-eq-npt-10ns.in -o equil-npt.out -p design_2.parm7 -c md_heat-300K.rst7 -r md_eq-npt.rst7 -x equil-npt.nc -ref md_heat-300K.rst7
