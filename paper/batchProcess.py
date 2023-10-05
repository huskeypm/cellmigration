#!/usr/bin/env python
import matplotlib.pylab as plt
import numpy as np
import pytraj as pt
import yaml 
import pandas as pd


# In[5]:


# insert and load the names of the pdb/dcd files that you made using the yamlFile to change the parameters
# calculate the mean square displacements for each pkl file
# plot the simulated trajectories of the particles
'''
Probably want to merge in processYaml to grab job info 

have try run meant to just get names of jobs, then do processing 
'''


cases=dict()
class empty:pass

#######
# user defined for now 
path="/home/pkh-lab-shared/migration/231004/"

yamlNamesFig4 =[
        "out_nParticles10.000000_00.yaml", #
"out_nParticles10.000000_01.yaml", #
"out_nParticles10.000000_02.yaml", #
"out_nParticles20.000000_00.yaml", #
"out_nParticles20.000000_01.yaml", #
"out_nParticles20.000000_02.yaml", #
"out_nParticles40.000000_00.yaml", #
"out_nParticles40.000000_01.yaml", #
"out_nParticles40.000000_02.yaml", #
"out_cellRad0.050000_00.yaml", #
"out_cellRad0.050000_01.yaml", #
"out_cellRad0.050000_02.yaml", #
"out_cellRad0.100000_00.yaml", #
"out_cellRad0.100000_01.yaml", #
"out_cellRad0.100000_02.yaml", #
"out_cellRad0.200000_00.yaml", #
"out_cellRad0.200000_01.yaml", #
"out_cellRad0.200000_02.yaml", #
"out_cellAttr0.500000_00.yaml", #
"out_cellAttr0.500000_01.yaml", #
"out_cellAttr0.500000_02.yaml", #
"out_cellAttr1.000000_00.yaml", #
"out_cellAttr1.000000_01.yaml", #
"out_cellAttr1.000000_02.yaml", #
"out_cellAttr2.000000_00.yaml", #
"out_cellAttr2.000000_01.yaml", #
"out_cellAttr2.000000_02.yaml"
]

yamlNamesFig5=[
"out_nCrowders8.000000_00.yaml",     
"out_nCrowders8.000000_01.yaml",     
"out_nCrowders8.000000_02.yaml",     
"out_nCrowders16.000000_00.yaml",     
"out_nCrowders16.000000_01.yaml",     
"out_nCrowders16.000000_02.yaml",     
"out_nCrowders32.000000_00.yaml",     
"out_nCrowders32.000000_01.yaml",     
"out_nCrowders32.000000_02.yaml",     
"out_crowderRad5.000000_00.yaml",     
"out_crowderRad5.000000_01.yaml",     
"out_crowderRad5.000000_02.yaml",     
"out_crowderRad10.000000_00.yaml",     
"out_crowderRad10.000000_01.yaml",     
"out_crowderRad10.000000_02.yaml",     
"out_crowderRad20.000000_00.yaml",     
"out_crowderRad20.000000_01.yaml",     
"out_crowderRad20.000000_02.yaml",     
        ]

########
equilFrame = 400
equilFrame = 0
dt = 1.

## 
## FUNC
##

# casename should include full path 
def ProcessTraj(caseName): 
    # load
    try:
      dcd=caseName+".dcd"; pdb=caseName+".pdb"
      traj = pt.iterload(dcd,pdb)
    except:
      raise RuntimeError("You're likely missing a file like %s"%dcd)

    ## get flux
    # for each particle, get dx in all directions, provide dt as input
    # select particles in some neighborhood of y=0?
    # grab dx along flux direction
    # sum(dx) / delta y

    ## get rmsd
    mask='@RC'
    rmsdAll = pt.rmsd(traj, mask='@RC', ref=0)
    rmsd = rmsdAll[equilFrame:]
    tEnd = np.shape(rmsd)[0]
    ts = np.arange(tEnd) * dt

    # fit for D
    #slope,intercept= np.polyfit(ts, ts*iter, 1)
    slope,intercept= np.polyfit(ts, rmsd, 1)
    #print(slope)
    #plt.plot(rmsd)
    #plt.plot(ts,ts*slope+intercept)
    #

    # x**2 = 4*D*t
    #      = slope * t ==> D = slope/4.
    Di = slope/4.
    print(Di) 
    return Di 

##
## MAIN 
## 
def doit(figName,yamlNames): 
#trajNames=[]
#tags=[]
  
  df = pd.DataFrame(columns=["trajName","tag","condVal","D"]) 
  for yamlName in yamlNames:
      # open yaml
      yamlName = path+"/"+yamlName
      with open(yamlName, 'r') as file:                                
        auxParams = yaml.safe_load(file)
      # get output name 
      trajName = auxParams['outName']                   
      print(trajName) 
      #trajNames.append(trajName) 
      # get 'tag'
      tag = auxParams['tag']
      # get cond value 
      condVal = auxParams[tag]
      print(tag,condVal)
  
      # process 
      Di = ProcessTraj(trajName) 
  
      # add to dataframe 
      df.loc[len(df.index)] = [trajName, tag, condVal,Di]
  
  print(df)
  
  outCsv = path+figName+".csv"       
  print("printed %s"%outCsv)
  df.to_csv(outCsv)               
  

import sys


#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -validation" % (scriptName)
  msg+="""
  
 
Notes:

"""
  return msg

#
# MAIN routine executed when launching this script from command line 
#
if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  #fileIn= sys.argv[1]
  #if(len(sys.argv)==3):
  #  1
  #  #print "arg"

  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-fig4"):
      figName = arg.replace("-","")
      yamlNames = yamlNamesFig4
      doit(figName,yamlNames) 
      quit() 

    elif(arg=="-fig5"):
      figName = arg.replace("-","")
      yamlNames = yamlNamesFig5
      doit(figName,yamlNames) 
      quit()

  





  raise RuntimeError("Arguments not understood")




