THis code was from a previous version. Deprecated

##############################################
## Other examples (Deprecated for now) 

WARNING: operator.py is now in Archive. I was getting a strange error about tuples that I nartrowed down to this script. 

Test case: 
```
$ python3 operator.py
```
or
```
$ /usr/bin/python3 operator.py
```

For the lengthy simulation, 

```
$ nohup python3 operator.py
```

Cases for final data (see notes below for additional instructions) 
```
# mutiple channel witdths 
$ python3 operatorFinal.py -pathwidths 

# open box sims 
```
## Configuration Guide for Various Cases
### Case 1: Open vs. Channel Model 
In this script, there are two ways to simulate the migration. 
1. To assess the motility (random motion) of cells, so called "Open" model is mimicking the cells placed in the large pool. To activate "open" model, the variable called "simType" needs to be set as "box". In ```operator.py```, this variable can be reassigned by the line with "Type" variable. As a default, it is "box" to simulate "open" model. If "Type" is given in "operator.py", there are only two options: "'box'" or "'non-box'". Any command other than "box" will activate the simulation with "slab", which is equivalent to "channel" structure (described in "runner.py".) 
2. Size Control: "UnitLength" represents the size of "box" or "channel" model. For instance, if you set "UnitLength = 1" then in the box model, it will create a square box with length of 10. "channel" model is a bit complicated. Due to the placement of channel between two reserviors or large pools, it takes another variable called "Indentation", which will control the width of "channel". Therefore, for the simplistic approach, it can be handy to adjust "UnitLength" but the downside is that it will result in differentiating the number of particles (cells) in the system. The best approach is to fix "UnitLength" and adjust "Indentation" with given allowance. For instance, "UnitLength = 50" then the reservior (square shape) length is "10". The channel width is determined by "Indentation", which is indicating how low from the height of reservior box and how high from the bottom of reservior box. If "Indentation = 1" then the width of channel is 8. For this reason, it is important "UnitLength" is substantially large enough to accomodate "Indentation" in the case where the simulation varies the width of channel. 
### Case 2: Different Densities
Densities are directly controlled by the variables called "Density1" and "Density2". 
1. Density1: density of cells in the reservior 1 (channel and box). This represent the number of "resting" or "unstimulated" cells. 
2. Density2: density of cells in the channel (not applicable for Box model). It is either all "resting" or "activated" 
3. numOfDeadCell: determins a number of crowder particles in the channel domain. It is also not applicable for the box model. The particle is larger than regular cells and they are placed with fixed distance from their neighboring crowder particle. 


## Analysis Guide 
### Generating CSV file for the RMSD/MSD analysis
The analysis is done by VMD script 
To run vmd without display or GUI, 
execute the following command

```
$ vmd -dispdev 
```

within the vmd environment, 

