#!/usr/bin/env python
import matplotlib.pylab as plt
import numpy as np

# pip3 install pytraj
import pytraj as pt


# In[5]:


# insert and load the names of the pkl files that you made using the yamlFile to change the parameters
# calculate the mean square displacements for each pkl file
# plot the simulated trajectories of the particles


cases=dict()
class empty:pass

#######
# user defined for now 
path="/home/pkh-lab-shared/migration/231004/"
condVals =[10,20,40]
trajNames = [
        "expt_nParticles10.000000",
        "expt_nParticles20.000000",
        "expt_nParticles40.000000"
        ]
iters=3
########
equilFrame = 400
equilFrame = 0
dt = 1.

import pandas as pd


slopes = []
Ds= []
Dstds= []
for trajName in trajNames:
  for iter in range(iters):

    caseName = path+trajName + "_%.2d"%(iter)
    print(caseName)
  
    # load
    try:
      dcd=caseName+".dcd"; pdb=caseName+".pdb"
      traj = pt.iterload(dcd,pdb)
    except:
      raise RuntimeError("You're likely missing a file like %s"%dcd) 

    # get rmsd 
    mask='@RC'
    rmsdAll = pt.rmsd(traj, mask='@RC', ref=0)
    rmsd = rmsdAll[equilFrame:]
    tEnd = np.shape(rmsd)[0] 
    ts = np.arange(tEnd) * dt
    
    # fit for D 
    #slope,intercept= np.polyfit(ts, ts*iter, 1)
    slope,intercept= np.polyfit(ts, rmsd, 1)
    slopes.append(slope)
    #print(slope)
    #plt.plot(rmsd)
    #plt.plot(ts,ts*slope+intercept)
    # 

  # x**2 = 4*D*t
  #      = slope * t ==> D = slope/4.
  iDs =  np.asarray( slopes )/4.

  case = empty()
  case.Dmean = np.mean( iDs )
  case.Dstd = np.std( iDs)
  Ds.append(case.Dmean)
  Dstds.append(case.Dstd)         

  print("Might phase out cases") 
  cases[trajName] = case

df = pd.DataFrame(
        {'trajName':trajNames,
         'nParticles':condVals,  
         'D':Ds,
         'Dstd':Dstds
        })

outCsv = path+"/test.csv"
print("printed %s"%outCsv)
df.to_csv(outCsv)               


