# cellmigration
1. [ Description ](#desc)
2. [ Getting started ](#start)
   [ Installation ](#install)
3. [ Usage tips ](#usage)
4. [ Analysis ](#analysis)
5. [ Advanced ](#advanced)


<a name="desc"></a>
# Description
- An openmm-based cell simulator
- Originally based on the OpenmmKant repo created by Ben Chun
- Code for running langevin particle simulations

# Getting started
<a name="start"></a>
- Refer to https://colab.research.google.com/drive/1Txu8ufScId-AqEvYnLM-yJx1kzMM_256?usp=sharing for a colab notebook that can quickly get you started
- The notebook is also copied here: conda_openmm_pytraj.ipynb


<a name="install"></a>
# Installation
## Python packages via anaconda
- See http://docs.openmm.org/latest/userguide/application/01_getting_started.html
- Create a new environment (seems cleaner) 
```
conda create -n openmm-env -c conda-forge openmm pyyaml matplotlib 
conda install -c conda-forge ambertools=23  
```

## from CLI 
- Check out the code 
```
git clone git@github.com:huskeypm/cellmigration.git
```

## To initialize environment 

if using conda 
- Activate env
```
conda activate openmm-env
```
- Test the installation 
```
python3 -c "import openmm"
```
- check to make sure nvidia is being used (if running on faust) 
```
nvidia-smi
```

if using kant
- Set env variable
```
source config.bash 
```

<a name="usage"></a>
## Execution 
- It is recommended to run the brownian .py from the command line via 
```
python3 brownian_v3.py -validation 
```

- To see the list of parameters used:
```
python3 brownian_v3.py -printVar
```


- The program is customized using parameters that are loaded in 'yaml' format. The syntax for calling the code with FILE.yaml is
```
python3 brownian_v3.py -yamlFile FILE.yaml -run
```
- Example yaml files are provided in the tests/ directory. Check out the README.md file in that directory for more details. 
- Run files and README.md for HPC (faust) are available in the ./run subdirectory

<a name="analysis"></a>
## Analysis
- Trajectory files like test.pkl can be opened and analyzed using the notebook bd_sims.ipynb in ./tests. Note that an example for computing mean square displacements (MSD) is provided therein. 
- code will also generate pdb/dcd files. These can be opened using vmd
-- VMD: load test.pdb. Right click entry in GUI. Select Load Data Into Molecule. Select dcd
- Lastly, if you know what you're doing, you can use the processYaml file
```
python processYaml.py -yaml tests/expt.yaml
```
This will print all of the parameters in csv format as well as an output yaml file


## Advanced
<a name="advanced"></a>

### Visualization
- Recommend using vmd to visusalize the output pdb and dcd files
  
### Yaml files
- Refer to tests/README.md for example inputs
  
### jupyter SSH 
- on remote host 
```
jupyter notebook --no-browser # --port 8888
```
- on local host 
```
ssh -L localhost:8890:localhost:8888    pkekeneshuskey@kant.luc.edu
```

### Fitting procedure (need to update) 
I adjusted the nUpdates parameter to equal the number of frames taken by the microscope
The framerate parameter is set to #/min 1 fr/90s  
The distance units in the code are assumed to be [um] though openmm assumes [nm]
The friction parameter was adjusted s.t. the MSD at the last frame was close to the expt value

### generating jobs for faust
Use master.ipynb on faust 

## TODO
- Can't get particle packing greater than 9/check that effectiveDim is helpful 
- ATP gradient is small
- Can we restrict particles to be in boundary left of the box 
- fermi dirac potential from Y 
- RESOLVED There is some kind of problem with the z-constraint. Compare brown_wnonbond.py z coordinates relative to openmuller.py in the tests directory 
- CANCEL PBCs in Y direction should be added 
- DONE implement flux calculation 
- DONE piecewise continuous functions? (this doesn't appear to be supported 
- DONE program fails with crowderDim>100



## Previous 
- Old notes are contained in README_v1.md
