Workflow
- Determine your code/data locations on the remote host (faust)
- in faust.bash, set run directory, runfile, and exec
 export SRC=/home/pkekeneshuskey/source/cellmigration/
 export PATH=$PATH:$SRC
 export EXEC="$PYTHON $SRC/brownian_v3.py"
 export RUNDIR=/data/pkekeneshuskey/231110/
 export RUNFILE="testrun"  # in RUNDIR  

- Do the same for local host 
 export LOCALDIR=/home/pkh-lab-shared/migration/231110/


- Generate jobs on kafka/kant (for now) using master.ipynb.  Be sure to set the outdir and cmd in the beginning of the notebook. Lines will be printed out to screen that should be pasted in a file like 'run' in the LOCALDIR
 $EXEC nocrowder_nCells10.000000_00.yaml -run
 $EXEC nocrowder_nCells10.000000_01.yaml -run
 $EXEC nocrowder_nCells10.000000_02.yaml -run

-rscp data to to remote host (faust). These are stored in a directory eg 
 scp -r $LOCALDIR/* $FAUST:$RUNDIR/                   # */     

- be sure that the appropriate RUNFILE name is used in the faust submission file 
! Note: suggested to split the runfiles into subsets so can run on multiple nodes (i use split_files.pl)

- jobs can be processed using batchProcess.py: 
python /home/pkekeneshuskey/source/cellmigration/batchProcess.py  -all



