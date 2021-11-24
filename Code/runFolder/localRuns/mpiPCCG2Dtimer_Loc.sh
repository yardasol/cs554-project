#!/bin/bash

n=11	## Size of 2D Poisson matrix = (n^2)x(n^2) 
p=5 	## No. of MPI ranks to use
nf=10	## no. of CG Full solves done to avg time
nc=10	## no. of CG CSR solves done to avg time
nd=10	## no. of CG DIA solves done to avg time

mpicc ../../Drivers/mpiPCCG2Dtimer.c -lm -o 2D_pccg_timer.exe
echo $'mpiCG2Dtimer.c Compiled with mpicc \n'

mpirun --oversubscribe -np $p 2D_pccg_timer.exe $n $nf $nc $nd
echo $'Job Run using mpirun - completed \n'

rm 2D_pccg_timer.exe
