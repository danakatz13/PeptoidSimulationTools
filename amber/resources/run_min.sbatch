#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=2:00:00

module load amber/openmpi/intel/22.00

sander -O -i min.in -p design_1.parm7 -c design_1.rst7 -o min1.out -r min1.rst
