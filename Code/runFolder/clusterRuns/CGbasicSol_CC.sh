#!/bin/bash

module load openmpi
echo "MPI module loaded"

mpicc ../../Drivers/CGbasicSolver.c -lm -o cgbasic_sol.exe
echo "CGbasicSolver.c Compiled with mpicc"

sbatch ./CGbasSol_CC.sbatch 
echo "Job submitted via sbatch"

