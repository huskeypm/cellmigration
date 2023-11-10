# cellmigration
- An openmm-based cell simulator
- ''Originally based on the OpenmmKant repo created by Ben Chun''
- ''Code for running langevin particle simulations''


# Installation
## Python packages via anaconda
See http://docs.openmm.org/latest/userguide/application/01_getting_started.html
```
conda install -c conda-forge openmm
```

## from CLI 
- Check out the code 
```
git clone https://github.com/bending456/OpenMMKant
```

- Revise the file config.bash to include the following environmental variables

```
export GOPATH=${HOME}/go
export PATH=/usr/local/go/bin:${PATH}:${GOPATH}:/usr/local/go:${HOME}/anaconda3/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/bchun/OpenMM/lib:/home/bchun/OpenMMpy2/lib:/home/bchun/OpenMM:/home/bchun/OpenMMpy2
export JUPYTER_PATH=/home/bchun/anaconda3/bin
export PATH=/usr/local/cuda-10.1/bin:/usr/local/cuda-10.1/bin:/home/bchun/.local/bin:/usr/local/go/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/home/bchun/go:/usr/local/go:/home/bchun/anaconda3/bin:/home/bchun/.openmpi/bin:/home/bchun/.openmpi/bin
export LD_LIBRARY_PATH=/usr/local/cuda-10.1/lib64:/usr/local/cuda-10.1/lib64::/home/bchun/OpenMM/lib:/home/bchun/OpenMMpy2/lib:/home/bchun/OpenMM:/home/bchun/OpenMMpy2:/home/bchun/.openmpi/lib/:/home/bchun/.openmpi/lib/
```

- add config.bash to your environment
```
source config.bash
```

- Test the installation 
```
python3 -c "import simtk"
```

## Execution 
- It is recommended to run brown_wnonbond.py from the command line via 
```
python3 brown_wnonbond.py -validation 
```

- To see the list of parameters used:
```
python3 brown_wnonbond.py -printVar
```


- The program is customized using parameters that are loaded in 'yaml' format. The syntax for calling the code with FILE.yaml is
```
python3 brown_wnonbond.py -yamlFile FILE.yaml -run
```

- An example yaml file is provided [here](https://github.com/bending456/OpenMMKant/blob/main/tests/paramSet1.yaml). In this example, the trajectory file is written to x.pkl


- Note: some of our installations are old, so you make have to import simtk.openmm instead of just openmm. If so, edit tests/brown_wnonbond.py accordingly


## Analysis
- Trajectory files like test.pkl can be opened and analyzed using the notebook bd_sims.ipynb in ./tests. Note that an example for computing mean square displacements (MSD) is provided therein. 
- code will also generate pdb/dcd files. These can be opened using vmd
-- VMD: load test.pdb. Right click entry in GUI. Select Load Data Into Molecule. Select dcd
- Lastly, if you know what you're doing, you can use the processYaml file
```
 python processYaml.py -yaml tests/expt.yaml
```
This will print all of the parameters in csv format as well as an output yaml file


## jupyter SSH 
  jupyter notebook --no-browser # --port 8888
ssh -L localhost:8890:localhost:8888    pkekeneshuskey@kant.luc.edu

## TODO
- RESOLVED There is some kind of problem with the z-constraint. Compare brown_wnonbond.py z coordinates relative to openmuller.py in the tests directory 
- CANCEL PBCs in Y direction should be added 
- DONE implement flux calculation 
- DONE piecewise continuous functions? (this doesn't appear to be supported 
- DONE program fails with crowderDim>100


## Fitting procedure (need to update) 
I adjusted the nUpdates parameter to equal the number of frames taken by the microscope
The framerate parameter is set to #/min 1 fr/90s  
The distance units in the code are assumed to be [um] though openmm assumes [nm]
The friction parameter was adjusted s.t. the MSD at the last frame was close to the expt value


Old notes are contained in README_v1.md
