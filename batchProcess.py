#!/usr/bin/env python
import matplotlib.pylab as plt
import numpy as np
import brown_util as bu

# pip3 install pytraj
import pytraj as pt


# In[5]:


# insert and load the names of the pkl files that you made using the yamlFile to change the parameters
# calculate the mean square displacements for each pkl file
# plot the simulated trajectories of the particles


cases=dict()

class empty:pass

trajNames = ["10","20"]
trajNames=["10"]
#use dcd
equilFrame = 400
equilFrame = 0
iters=1

slopes = []
for trajName in trajNames:
  for iter in range(iters):

    caseName = trajName + "_%d"%(iter)
    print(caseName)
  
    # load
    try:
      traj = pt.iterload(caseName+".dcd", caseName+".pdb")
    except:
      RaiseRuntimeError("You're likely missing a file ") 

    # get rmsd 
    mask='@RC'
    rmsdAll = pt.rmsd(traj, mask='@RC', ref=0)
    rmsd = rmsdAll[equilFrame:]
    tEnd = np.shape(rmsd)[0] 
    dt = 1.
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
  Ds =  np.asarray( slopes )/4.

  case = empty()
  case.Dmean = np.mean( Ds )
  case.Dstd = np.std( Ds)

  cases[trajName] = case
quit()


