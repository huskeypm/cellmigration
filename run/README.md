Workflow

Generating jobs on kafka/kant and rscp to faust 
- these are stored in a directory eg 
scp -r migration_data/231110/ $FAUST:/data/pkekeneshuskey/

- copy of run files 
scp -r run* $FAUST:/home/pkekeneshuskey/source/cellmigration/

- adjust paths in run file if needed
---> /data/pkekeneshuskey/231110/XXXX.yaml

- in faust.bash
replace RUNFILE with appropriate name 
--> run_fig4

!!!! need to change paths in yamlFiles to refer to local directory only 
same thing for run files 





