#!/bin/bash

n=11	## Size of 1D Poisson matrix
p=5 	## No. of MPI ranks to use
nf=10	## no. of CG Full solves done to avg time
nc=10	## no. of CG CSR solves done to avg time
nd=10	## no. of CG DIA solves done to avg time

mpicc ../../Drivers/mpiCG1Dtimer.c -lm -o 1Dtimer.exe
echo $'mpiCG1Dtimer.c Compiled with mpicc \n'

mpirun --oversubscribe -np $p 1Dtimer.exe $n $nf $nc $nd
echo $'Job Run using mpirun - completed \n'

rm 1Dtimer.exe
