#!/bin/bash

module load openmpi
echo "MPI module loaded"

mpicc ../../Drivers/mpiCG2Dtimer.c -lm -o 2Dtimer.exe
echo "mpiCG2Dtimer.c Compiled with mpicc"

sbatch ./mpiCG2Dtim_CC.sbatch 
echo "Job submitted via sbatch"
