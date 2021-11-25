#!/bin/bash

n=11	## Size of 1D Poisson matrix
p=5 	## No. of MPI ranks to use

mpicc ../../Drivers/PCCGbasicSolver.c -lm -o pccgbasic_sol.exe 
echo $'PCCGbasicSolver.c Compiled with mpicc \n'

mpirun --oversubscribe -np $p pccgbasic_sol.exe $n
echo $'Job Run using mpirun - completed \n'

rm pccgbasic_sol.exe
