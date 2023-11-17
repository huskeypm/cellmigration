#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5 
#SBATCH --gres=gpu:1
#SBATCH --time=1:00:00
#SBATCH --job-name=cam1Equi
#SBATCH --mem=40G

# load modules 
module load shared
module load cuda11.3/toolkit/11.3.0

# this was a conda instance that includes openmm
conda activate SEEKR

#export PATH="$PATH:/home/pkekeneshuskey/source/cellmigration/"


# assumes paths, etc are correct 
# See readme 
# exec  ~/source/cellmigration/brownian_v3.py
RUNDIR=/data/pkekeneshuskey/231110/                     
RUNFILE="testrun"  # in RUNDIR  
cd $RUNDIR 
#bash $RUNFILE
#python -c "import openmm; print('imported')" 
pwd 
bash $RUNFILE
#python3 /home/pkekeneshuskey/source/cellmigration/brownian_v3.py 

# WARNING: need to copy files to TEMP/copy BACK

echo "FINISHED" 

