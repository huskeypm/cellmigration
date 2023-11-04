Use writeInputs for creating appropriate yaml files from templates 

* Setup
make a reference file that includes ATP and crowders
 templatefull.yaml

* Batch processing 
edit master.ipynb to run/process jobs
upload/process files with colab notebook above 




Instructions are inside the script 

Simulation Procedure 
* Run writeInputs.py without inputs
 python3 writeInputs.py -fig[4,5,???] | grep python > runallXXXX
It will generate output yaml files

* Run each using this setup: 
 source ../config.bash
 python3 brown_wnonbond.py -yamlFile paper/out_nParticles10.000000.yaml -run

Note: this is spit out by the program anyway, so pipe to a 'runall' bash script, so that you can run
 bash runallXXXX

Analysis Procedure:
* Need to edit the inputs to this to make sure data is processed correctly. Namely, need to define the yamlNamesFig4 var at the beginning of the script
 source ../config.bash
 python3 batchProcess.py -fig[4,5,??]

* using colab for now for analysis
Upload <case>.csv from output directory to google colab
https://colab.research.google.com/drive/1AZdQw3q9cdjF7HGfCZCfOkKR7idzIYLD#scrollTo=8uJckvVigVb9

(if changes are made to the colab doc, be sure to download and commit to repo) 



exptl.yaml <- reproduces experimentally-observed diffusion
