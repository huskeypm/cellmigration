Yaml files are used to feed parameters into the simulation. A basic example is provided below:
```
# example.yaml
nCells: 30  # num of cells 
cellRad:  1   # cell radii 
cellAttr: 0.0   # attraction constant 
nCrowders:  1   # number of crowedrs
crowderRad: 1   # um
crowderAttr: 0.0

domainXDim: 200  # size of simulation domain
domainYDim: 200  # size of simulation domain
nUpdates: 500 # number of steps            
friction:  10 # [1/ps] friction coefficient 

containmentPotential: None     #  to contain particles 
xPotential: False   # [True/False]  chemotraction 
xScale: 0.1        # scale for chemoattractant gradient (when xpotential is true)  
yPotential: False # [True/False] chemotraction
outName: expt
```
Examples in this directory
- square.yaml - imposes square (rectangle actually) potential to contain particles
- cylinder.yaml - imposes cylinder potential to contain particles
- expt.yaml - parameters needed to reproduce expt w/o ATP 
- expt_watp.yaml - parameters needed to reproduce expt w ATP 
- crowder_gradient.yaml - finite sized crowder with chemokine gradients in two directions
- effrad.yaml - assumes crowders have an effective radius based on running a probability calculation 
- states.yaml - instructions for simulating microglia states 
