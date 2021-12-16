#!/bin/bash

module load openmpi
echo "MPI module loaded"
TASKS=(1 2 4 8 16 32 64)
PSIZE=(128 256 512 1024 2048)
TP=("1D")
FPRE="PCCG"
FPOST="timer"
for tp in "${TP[@]}"
do
    DIR="./$tp"
    cd $DIR
    mpicc ../../../Drivers/mpiPCCG$tp$FPOST.c -lm -o $tp$FPRE$FPOST.exe
    echo "mpiPCCG$tp$FPOST.c Compiled with mpicc"
    for s in "${PSIZE[@]}"
    do
	for t in "${TASKS[@]}" 
	do

            F="s$s-t$t.sbatch"
	    sbatch $F
	    echo "job submitted via sbatch"
        done
    done
    cd ../
done
