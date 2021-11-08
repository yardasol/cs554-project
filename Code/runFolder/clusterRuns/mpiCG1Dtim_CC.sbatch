#!/bin/bash

#####################################################
##        Important parameters of your job         ##
##             are specified here 		   ##
#####################################################

#SBATCH -J CG_1D_Timer				## job name
#SBATCH -t 02:50:00					## total computing time
#SBATCH -N 2						## number of nodes 
#SBATCH -n 33						## total number of MPI tasks
####SBATCH --mem-per-cpu=7200			## total memory required (if not specified, will get the entire node memory)
#SBATCH --partition=cs				## queue
#SBATCH --output=../Results/1Dtimer.out%j		## file that will receive output from execution
#SBATCH --error=../Results/1Dtimer.err%j			## file that will receive any error messages

#SBATCH --mail-user=gautham3@illinois.edu
#SBATCH --mail-type=begin 			## email me when the job starts
#SBATCH --mail-type=end 			## email me when the job finishes


########## Run your executable ######################
n=8192			## Size of 1D Poisson matrix
nf=1			## no. of CG Full solves done to avg time
nc=1			## no. of CG CSR solves done to avg time
nd=1			## no. of CG DIA solves done to avg time

module load openmpi
srun ./1Dtimer.exe $n $nf $nc $nd
