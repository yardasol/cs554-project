#!/bin/bash

module load openmpi
echo "MPI module loaded"

mpicc ../../Drivers/mpiPCCG1Dtimer.c -lm -o 1DPCCGtimer.exe
echo "mpiPCCG1Dtimer.c Compiled with mpicc"

sbatch ./mpiPCCG1Dtim_CC.sbatch 
echo "Job submitted via sbatch"
