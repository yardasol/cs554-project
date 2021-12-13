
TASKS=(1 2 4 8 16 32 64)
PSIZE=(128 256 512 1024 2048)
NODES=4
TP=("1D" "2D")
TPNAME="PCCGtimer"
for tp in "${TP[@]}"
do
    DIR="./$tp"
    for s in "${PSIZE[@]}"
    do
	for t in "${TASKS[@]}" 
	do
	    if [ $t -lt $NODES ]
	    then
	        N=$t
	    else
	        N=$NODES
	    fi
	    echo "s$s-t$t.sbatch"
	    F=$DIR/s$s-t$t.sbatch
	    touch $F
	    echo "#!/bin/bash" | cat - > $F
	    echo "#SBATCH -J ILU-$tp-Timers-$s-t$t" | cat - >> $F
	    echo "#SBATCH --time=01:50:00" | cat - >> $F
	    echo "#SBATCH -N $N" | cat - >> $F
	    echo "#SBATCH -n $t" | cat - >> $F
	    echo "#SBATCH --partition=cs" | cat - >> $F
	    echo "#SBATCH --output=../../Results/$tp$TPNAME-s$s-t$t.out" | cat - >> $F
	    echo "#SBATCH --error=../../Results/$tp$TPNAME-s$s-t$t.err" | cat - >> $F
	    echo "#SBATCH --mail-user=oyardas2@illinois.edu" | cat - >> $F
	    echo "#SBATCH --mail-type=begin" | cat - >> $F
	    echo "#SBATCH --mail-type=end" | cat - >> $F

	    echo "n=$s" | cat - >> $F
	    echo "nf=100" | cat - >> $F
	    echo "nc=1000" | cat - >> $F
	    echo "nd=1000" | cat - >> $F
	    echo "module load openmpi" | cat - >> $F
	    echo "srun ./$tp$TPNAME.exe \$n \$nf \$nc \$nd" | cat - >> $F
        done
    done
done
