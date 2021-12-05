#!/bin/bash

module load openmpi
echo "MPI module loaded"

mpicc ../../Drivers/mpiPCCG2Dtimer.c -lm -o 2DPCCGtimer.exe
echo "mpiPCCG2Dtimer.c Compiled with mpicc"

sbatch ./mpiPCCG2Dtim_CC.sbatch 
echo "Job submitted via sbatch"
