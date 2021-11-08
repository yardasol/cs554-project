#!/bin/bash

n=11	## Size of 1D Poisson matrix
p=5 	## No. of MPI ranks to use

mpicc ../../Drivers/CGbasicSolver.c -lm -o cgbasic_sol.exe 
echo $'CGbasicSolver.c Compiled with mpicc \n'

mpirun --oversubscribe -np $p cgbasic_sol.exe $n
echo $'Job Run using mpirun - completed \n'

rm cgbasic_sol.exe
