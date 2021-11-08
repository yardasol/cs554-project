#!/bin/bash


n=11	## Size of 1D Poisson matrix

#####################################################
##        Important parameters of your job         ##
##             are specified here 		   ##
#####################################################

#SBATCH -J CG_1Dbasic_Sol			## job name
#SBATCH -t 00:50:00					## total computing time
#SBATCH -N 1						## number of nodes 
#SBATCH -n 5						## number of tasks per node
####SBATCH --mem-per-cpu=7200			## total memory required (if not specified, will get the entire node memory)
#SBATCH --partition=cs				## queue
#SBATCH --output=../Results/cgbsol.out%j		## file that will receive output from execution
#SBATCH --error=../Results/cgbsol.err%j			## file that will receive any error messages

#SBATCH --mail-user=gautham3@illinois.edu
#SBATCH --mail-type=begin 			## email me when the job starts
#SBATCH --mail-type=end 			## email me when the job finishes


########## Run your executable ######################
module load openmpi
srun ./cgbasic_sol.exe $n 
