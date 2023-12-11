WARNING: for the time being, use SEEKER conda, not openmm-env w pytraj
SEEKRNew



Workflow
- Determine your code/data locations on the remote host (faust)
- in faust.bash, set run directory, runfile, and exec
 export SRC=/home/pkekeneshuskey/source/cellmigration/
 export PATH=$PATH:$SRC
 export EXEC="$PYTHON $SRC/brownian_v4.py"   #  DOUBLE CHECK SINCE THIS GETS UPDATED FREQUENTLY!
 export RUNDIR=/data/pkekeneshuskey/231129/
 export RUNFILE="testrun"  # in RUNDIR  

- Do the same for local host (if generating jobs there)  
 export LOCALDIR=/home/pkekeneshuskey/data/231129/
 #export LOCALDIR=/home/pkh-lab-shared/migration/231129/


- Generate jobs on kafka/kant (or faust now) using master.ipynb.  Be sure to set the outdir and cmd in the beginning of the notebook. Lines will be printed out to screen that should be pasted in a file like 'run' in the LOCALDIR
 $EXEC nocrowder_nCells10.000000_00.yaml -run
 $EXEC nocrowder_nCells10.000000_01.yaml -run
 $EXEC nocrowder_nCells10.000000_02.yaml -run


or, now we can define a 'writeFile'; the contents of which can be cat to a runfile. In this example I wrote
 cat *sh > master   # */ 

- be sure that the appropriate RUNFILE name is used in the faust submission file 
! Note: suggested to split the runfiles (e.g. master from the prev example)  into subsets so can run on multiple nodes (i use split_files.pl)
 split_files.pl master 10
 cp *master $LOCALDIR/  #*


-rscp data to to remote host (faust). These are stored in a directory eg 
 ssh faust "mkdir $RUNDIR"
 rsync -rh $LOCALDIR/* $FAUST:$RUNDIR/                   # */     


- create run files and edit 
  cp faust_template.bash 0Xmaster.bash
  # edit 0Xmaster.bash to point to right file 

- execute
 sbatch faust_template.sh

- jobs can be processed using batchProcess.py (from dir where files were run): 
 python /home/pkekeneshuskey/source/cellmigration/batchProcess.py  -all

on faust, use 
  conda activate AmberTools23
  python -c 'import pytraj'

- after creating the output csv file in the previous step, the results can be processed via google colab
https://colab.research.google.com/drive/1AZdQw3q9cdjF7HGfCZCfOkKR7idzIYLD#scrollTo=TZf5qh69fbAw


.......
Launching a buncho jobs 

split_files.pl master 150
perl -pe 's/AA/01/g;' faust_template.bash > 01.bash



for i in {1..9}; do
    echo "$i"
    perl -pe "s/AA/0$i/g;" faust_template.bash > 0${i}_process.bash
    #sbatch 0${i}_process.bash
done
for i in {10..20}; do
    echo "$i"
    perl -pe "s/AA/$i/g;" faust_template.bash > ${i}_process.bash
    #sbatch ${i}_process.bash
done


find jobs that haven't been run 
-----
import os;
import glob
files = glob.glob('*.yaml')

for file in files:
  dcd=file.replace('yaml','dcd')
  if os.path.exists(dcd) is not True:
    print(file)
-------

rad > 12 isn't good
it seems like master is not doing the points I requested - double check 
crowder_atp_nCrowders9.000000_02.yaml <-- calling error about square rootable
too many cells 
figure out how to assemble all of the assembl;ed csv files 




