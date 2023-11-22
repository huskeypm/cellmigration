#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5 
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --job-name=openmmcellmigr
#SBATCH --mem=40G

# load modules 
module load shared
module load cuda11.3/toolkit/11.3.0

# this was a conda instance that includes openmm
conda activate SEEKR
export PYTHON=`which python` # /cm/shared/apps/anaconda3/envs/SEEKR/bin/python


# assumes paths, etc are correct 
# See readme 
export SRC=/home/pkekeneshuskey/source/cellmigration/
export EXEC="$PYTHON $SRC/brownian_v3.py"
export RUNDIR=/data/pkekeneshuskey/231117/                     
export RUNFILE="run_a"  # in RUNDIR  

cd $RUNDIR 
#python -c "import openmm; print('imported')" 
pwd 
bash $RUNFILE

# WARNING: need to copy files to TEMP/copy BACK

echo "FINISHED" 

