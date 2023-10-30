#!/usr/bin/env python
import matplotlib.pylab as plt
import numpy as np
import pytraj as pt
import yaml 
import pandas as pd
import glob


# In[5]:


# insert and load the names of the pdb/dcd files that you made using the yamlFile to change the parameters
# calculate the mean square displacements for each pkl file
# plot the simulated trajectories of the particles

cases=dict()
class empty:pass

#######
# user defined for now 
path="/home/pkh-lab-shared/migration/231004/"
#path=""

########
equilFrame = 400
equilFrame = 0
dt = 1.

## 
## FUNC
##
## get flux
def GetFlux(traj, mask='@RC',display=False):
# for each particle, get dx in all directions, provide dt as input
# select particles in some neighborhood of y=0?
# grab dx along flux direction
# sum(dx) / delta y
  numFrames = (np.shape(traj.xyz))[0]

  ## get cells 
  indices = pt.select_atoms(traj.top, mask) 
  
  # xyz: frames, natoms, coords
  #print(traj.xyz[2,:,0])
  xThresh = -0.

  
  plt.figure()
  diffs = traj.xyz[-1,indices,0] 
  diffs -= traj.xyz[0,indices,0] 
  #print(traj.xyz[0,indices,0])
  #print(traj.xyz[-1,indices,0])
  plt.hist(diffs)
  plt.gcf().savefig("diffs.png") 
  


  
  # below Thresh (only need to track one 'compartment' for now since we are dividing the space into two sections
  # along x direction 
  #l =np.array(traj.xyz[:,0,0] < xThresh, dtype=int) 
  #r =np.array(traj.xyz[:,0,0] >=xThresh, dtype=int) 
  #print(np.sum(l),np.sum(r))
  l =np.array(traj.xyz[:,indices,0] < xThresh, dtype=int) 
  l =np.sum(l,axis=1)
  # if the population in the compartment changes between two times, then a particle has entered/left
  diff = np.diff(l)  
  fluxArea = diff/dt # not normalizing by area
  JA = np.average(fluxArea)
  if (l[-1] < 0.1*l[0]):
      print("WARNING: compartment is nearly depleted/flux estimates may be unreliable") 
  #print("J*A = %f "%JA)
  
  #x = traj.xyz[:,indices,0]
  #y = traj.xyz[:,indices,1]
  #plt.plot(x,y)
  #plt.gcf().savefig("testxy.png") 
  if display:
    print(l[0],l[-1],np.sum(fluxArea),np.average(fluxArea)*numFrames) # 320 frames 
    plt.plot(l,label="#particles in x<thresh")     
    plt.plot(fluxArea,label="flux*area")
    plt.legend(loc=0)
    plt.gcf().savefig("test.png") 

  return JA         

def CalcD(traj,mask='@RC',csvName=None):
  rmsdAll = pt.rmsd(traj, mask='@RC', ref=0)
  rmsd = rmsdAll[equilFrame:]
  tEnd = np.shape(rmsd)[0]
  ts = np.arange(tEnd) * dt
  # fit for D
  slope,intercept= np.polyfit(ts, rmsd, 1)
  #print(slope)
  #plt.plot(rmsd)
  #plt.plot(ts,ts*slope+intercept)
  #

  if csvName is not None:
    dim = np.shape(rmsd)[0]
    csv = np.reshape(np.zeros(dim*2),[dim,2])
    csv[:,0] = ts
    csv[:,1] = rmsd
    np.savetxt(csvName+".csv",csv)   

    plt.figure()
    plt.plot(ts,rmsd)
    plt.gcf().savefig(csvName+".png")

  # x**2 = 4*D*t
  #      = slope * t ==> D = slope/4.
  Di = slope/4.
  return Di 

# assumes cylindrical inclusions 
def CalcVolFrac(auxParams): 
  d = auxParams['domainDim']
  n = auxParams['nCrowders']
  r = auxParams['crowderRad']

  volFrac = d*d - n*np.pi*r*r  # since domain is square
  print(d)
  print(n,r)
  print(n*np.pi*r*r)
  volFrac /= d*d                

  return volFrac 

# casename should include full path 
def ProcessTraj(caseName,display=False): 
    # load
    try:
      dcd=caseName+".dcd"; pdb=caseName+".pdb"
      traj = pt.iterload(dcd,pdb)
    except:
      raise RuntimeError("You're likely missing a file like %s"%dcd)
    print("Loaded %s"%dcd)

    ## get D    

    Di=CalcD(traj,mask='@RC',csvName=caseName)                        
    JA=GetFlux(traj,mask='@RC',display=display)
    print("Di %f J %f"%(Di,JA))
    return Di,JA 

##
## MAIN 
## 

# reads the default params and those in the yaml file 
def doit(figName,yamlNamePrefix="*",single=False): 

  # get names
  if single:
      yamlNames=[yamlNamePrefix+".yaml"]
  else: 
    globTag = path+"/"+yamlNamePrefix+'*yaml'
    try: 
      yamlNames = glob.glob(globTag)
    except:
      raise RuntimeError("No files found") 

  #print(globTag) 
  #print(yamlNames)
  
  df = pd.DataFrame(columns=["trajName","tag","condVal","D","flux*A"]) 

  for yamlName in yamlNames:
      # open yaml
      #yamlName = path+"/"+yamlName
      with open(yamlName, 'r') as file:                                
        auxParams = yaml.safe_load(file)
      # get output name 
      trajName = auxParams['outName']                   
      print(trajName) 

      # compute vol frac
      volFrac = CalcVolFrac(auxParams)
      print(volFrac)
      quit()

      #trajNames.append(trajName) 
      # get 'tag'
      # get cond value 
      try:
        tag = auxParams['tag']
        condVal = auxParams[tag]
      except:
        tag = 'tag'
        condVal=1.
      print(tag,condVal)
  
      # process 
      display = False
      if single:
          display=True
      Di,JA = ProcessTraj(trajName,display=display) 
  
      # add to dataframe 
      df.loc[len(df.index)] = [trajName, tag, condVal,Di,JA]
  
  if single:
      return Di,JA
  
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
  msg+="  %s -fig4/-fig5/-single [yaml/dcdprefix]" % (scriptName)
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
  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    if(arg=="-fig4"):
      figName = arg.replace("-","")
      doit(figName,yamlNamePrefix="noatp_")
      doit(figName+"atp",yamlNamePrefix="atp_") 
      quit() 

    elif(arg=="-fig5"):
      figName = arg.replace("-","")
      doit(figName,yamlNamePrefix="crwdnoatp_")   
      doit(figName+"atp",yamlNamePrefix="crwdatp_")
      quit()

    elif(arg=="-single"): 
      yamlName=sys.argv[i+1]#
      Di,JA = doit("test.png",yamlName,single =True)
      #Di,JA = ProcessTraj(trajName,display=True)
      quit()

  





  raise RuntimeError("Arguments not understood")




