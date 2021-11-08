#!/bin/bash

module load openmpi
echo "MPI module loaded"

mpicc ../../Drivers/mpiCG1Dtimer.c -lm -o 1Dtimer.exe
echo "mpiCG1Dtimer.c Compiled with mpicc"

sbatch ./mpiCG1Dtim_CC.sbatch 
echo "Job submitted via sbatch"
