#!/usr/bin/env python
import matplotlib.pylab as plt
import brown_util as bu
import numpy as np
import pytraj as pt
import yaml 
import pandas as pd
import glob

print("SHOULD PULL INTO BU") 


# In[5]:


# insert and load the names of the pdb/dcd files that you made using the yamlFile to change the parameters
# calculate the mean square displacements for each pkl file
# plot the simulated trajectories of the particles

cases=dict()
class empty:pass

#######
# user defined for now 
path="/home/pkh-lab-shared/migration/231004/"
path="/home/pkh-lab-shared/migration/231110/"
#path=""

########
equilFrame = 400
equilFrame = 0
FRAME_CONV = 0.1 # min/fr 
dt = FRAME_CONV  # [min] 

## 
## FUNC
##
def CalcRDF(traj,
            mask1, # solvent
            mask2, # solute
            space=0.1,
            bins=20):
    #radial = pt.rdf(traj, solvent_mask=':WAT@O', solute_mask=':WAT@O', bin_spacing=0.2, maximum=12.)
    radial = pt.rdf(traj, solvent_mask=mask1, solute_mask=mask2, bin_spacing=space, maximum=bins) 

    return radial 


def CalcProbDist(traj, mask='@RC',display=False):
# for each particle, get dx in all directions, provide dt as input
# select particles in some neighborhood of y=0?
# grab dx along flux direction
# sum(dx) / delta y
  numFrames = (np.shape(traj.xyz))[0]

  ## get cells
  indices = pt.select_atoms(traj.top, mask)

  xs = traj.xyz[0:,indices,0]
  xs = np.ndarray.flatten(xs)  # n particles x m timesteps 
  ys = traj.xyz[0:,indices,1]
  ys = np.ndarray.flatten(ys)  # n particles x m timesteps 

  p,x,y= np.histogram2d(xs,ys,density=True)
  X, Y = np.meshgrid(x, y)

  # get PMF
  kT = 1
  thresh=1e-9
  p[p<thresh]=thresh
  pmf = -np.log(p) * kT 
  pmf-=np.min(pmf) 


  #display=True
  if display:
    plt.figure()
    plt.axis('equal')
    plt.pcolormesh(X, Y, p.T) # probably T is appropriate here 
    plt.colorbar()
    plt.gcf().savefig("prob.png",dpi=300)

    plt.figure()
    plt.axis('equal')
    plt.pcolormesh(X, Y, pmf.T) # probably T is appropriate here 
    plt.colorbar()
    plt.gcf().savefig("pmf.png",dpi=300)

  return p,X,Y


## get flux
def CalcFlux(traj, mask='@RC',display=False):
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

  if display: 
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

def CalcD(traj,mask='@RC',csvName=None, display=False):
  # in A^2 
  rmsdAll = pt.rmsd(traj, mask='@RC', ref=0)
  rmsd = rmsdAll[equilFrame:]
  msd = np.sqrt(rmsd) 
  # in NM^2
  AA_to_NMNM=1e-2
  msd_NMNM = msd*AA_to_NMNM

  # in XXX min/fr
  tEnd = np.shape(rmsd)[0]
  ts_MIN = np.arange(tEnd) * dt

  # fit for D [nm^2/min]
  slope,intercept= np.polyfit(ts_MIN, msd_NMNM, 1)
  #print(slope)
  #plt.plot(rmsd)
  #plt.plot(ts,ts*slope+intercept)
  #

  if csvName is not None:
    dim = np.shape(rmsd)[0]
    csv = np.reshape(np.zeros(dim*2),[dim,2])
    csv[:,0] = ts_MIN 
    csv[:,1] = msd_NMNM
    np.savetxt(csvName+".csv",csv,header='time[min],msd[nmnm]')   

  if display: 
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

  volFrac = bu.CalcVolFrac(n,d,r)

  return volFrac 

# casename should include full path 
def LoadTraj(caseName,warningOnly=False):
    # load
    try:
      dcd=caseName+".dcd"; pdb=caseName+".pdb"
      traj = pt.iterload(dcd,pdb)
    except:
      if warningOnly:
        print("You're likely missing a file like %s"%dcd)
        return None 
      else: 
        raise RuntimeError("You're likely missing a file like %s"%dcd)
    print("Loaded %s"%dcd)
    return traj


def ProcessTraj(caseName,display=False): 
    # LoadTraj
    traj = LoadTraj(caseName)

    ## get J,D    
    Di=CalcD(traj,mask='@RC',csvName=caseName)                        
    JA=CalcFlux(traj,mask='@RC',display=display)
    print("Di %f J %f"%(Di,JA))

    # get 2D histogram of populations
    CalcProbDist(traj,mask='@RC',display=display)
    return Di,JA 

##
## MAIN 
## 

# reads the default params and those in the yaml file 
def processYamls(figName,
        path=path,
        yamlNamePrefix="*",
        prefixOptions=None,# can list cellAttr etc to fine tune search 
        single=False,display=True,warningOnly=False): 

  # get names
  if single:
      yamlNames=[yamlNamePrefix+".yaml"]
  else: 
    globTag = path+"/"+yamlNamePrefix+'*yaml'
    try: 
      yamlNames = glob.glob(globTag)
    except:
      raise RuntimeError("No files found") 

  if prefixOptions is not None:
      yamlNames=[]
      for opt in prefixOptions:
        globTag = path+"/"+yamlNamePrefix+'_%s*yaml'%opt
        yamlNames+=glob.glob(globTag)
        

  #print(globTag) 
  #print(yamlNames)
  
  df = pd.DataFrame(
          columns=["trajName","tag","condVal","D","flux*A","Vol Frac"]
          ) 

  #print(yamlNames[0]) 
  for yamlName in yamlNames:
      # open yaml
      #yamlName = path+"/"+yamlName
      with open(yamlName, 'r') as file:                                
        auxParams = yaml.safe_load(file)
      # get output name 
      trajName = auxParams['outName']                   
      #print(trajName) 

      # compute vol frac
      volFrac = CalcVolFrac(auxParams)

      #trajNames.append(trajName) 
      # get 'tag'
      # get cond value 
      try:
        tag = auxParams['tag']
        condVal = auxParams[tag]
      except:
        tag = 'tag'
        condVal=1.
      #print(tag,condVal)
  
      # process 
      if LoadTraj(trajName,warningOnly) is None:
        continue 
      
      Di,JA = ProcessTraj(trajName,display=display) 
  
      # add to dataframe 
      df.loc[len(df.index)] = [trajName, tag, condVal,Di,JA,volFrac]
  
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
  msg+="  %s -fig4/-fig5/-single [yaml/dcdprefix] -all" % (scriptName)
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
      processYamls(figName,yamlNamePrefix="noatp_")
      processYamls(figName+"atp",yamlNamePrefix="atp_") 
      quit() 

    elif(arg=="-fig5"):
      figName = arg.replace("-","")
      processYamls(figName,yamlNamePrefix="crwdnoatp_")   
      processYamls(figName+"atp",yamlNamePrefix="crwdatp_")
      quit()

    elif(arg=="-single"): 
      yamlName=sys.argv[i+1]#
      Di,JA = processYamls("test.png",yamlName,single =True,display=True)
      #Di,JA = ProcessTraj(trajName,display=True)
      quit()

    elif(arg=="-all"): 
      processYamls("test.png",display=False,warningOnly=True,path="./")
      quit()

  raise RuntimeError("Arguments not understood")




