#!/usr/bin/env python
import matplotlib.pylab as plt
import brown_util as bu
import numpy as np
import pytraj as pt
import yaml 
import pandas as pd
import glob

# In[5]:


# insert and load the names of the pdb/dcd files that you made using the yamlFile to change the parameters
# calculate the mean square displacements for each pkl file
# plot the simulated trajectories of the particles
##
## MAIN 
## 

cases=dict()
class empty:pass

#######
# user defined for now 
path="/home/pkh-lab-shared/migration/231004/"
path="/home/pkh-lab-shared/migration/231110/"
path="/home/pkh-lab-shared/migration/231117/"
#path=""

########
equilFrame = 400
equilFrame = 0
FRAME_CONV = 0.1 # min/fr 
dt = FRAME_CONV  # [min] 

crowderDim = 50 


def ProcessTraj(caseName,display=False,auxParams=None): 
 # LoadTraj
 caseName = caseName.replace('.yaml',"")
 traj = bu.LoadTraj(caseName)
 traj._force_load =  True

 #, get RDF
 bu.CalcRDFs(traj,"RC","AC")
    #mask1='@RC'
    #mask2=':21@AC'# last atom (crowder) 
    #bu.CalcRDF(traj,mask1=mask1,mask2=mask2)  
#
 ## get J,D    
 Di=bu.CalcD(traj,mask='@RC',csvName=caseName)                        
 #xThresh = 300 # AA
 #xThresh = 600 # AA
 #JA=bu.CalcFluxLine(traj,mask='@RC',display=display,xThresh=xThresh)
 #print("Di %f J %f"%(Di,JA))

 # get 2D histogram of populations
 resolution = 4  # [px/um ]
 if auxParams is not None:
   dimx = auxParams['domainXDim']*2 # since +/- domainXDim
   dimy = auxParams['domainYDim']*2 # since +/- domainXDim
   dimsPx = resolution * np.array([dimx,dimy])
   midx = int( dims[0]/2 ) 

   #reservoirGap = auxParams['domainXDim'] - auxParams['crowderXDim']
   crowderDimsPx=[midx - auxParams['crowderXDim']*resolution,
                   midx + auxParams['crowderXDim']*resolution]
   print("x=",dimsPx) 
   print(crowderDimsPx)

 else:
   dimsPx = None 
   crowderDimsPx = None 
   #print("DEBUGGING") 
   #dummy = bu.CalcProbDist(traj,mask='@RC',
   #            tMax=1000000, caseName=caseName+"_1M",
#		display=False,bins=dims)
#    dummy = bu.CalcProbDist(traj,mask='@RC',
#                tMax=3000000, caseName=caseName+"_3M",
#		display=False,bins=dims)
 prob,d,d,dx,dy = bu.CalcProbDist(traj,mask='@RC',caseName=caseName,display=display,bins=dimsPx)

 #print("WARNING: crowederDim is hard-coded") 
 #xLimCrowder=[500,1000]
 areaFrac,Jcrowd,Jreserv =  bu.CalcAverageFlux(
   prob,
   D=1,
   barrierLims = crowderDims,
   reservoirLims = [1000,1200],
   dx=dx,
   dy=dy )

 return Di,Jreserv, Jcrowd, areaFrac 

def GetParams(yamlName):
  with open(yamlName, 'r') as file:                                
    auxParams = yaml.safe_load(file)

  return auxParams 


# reads the default params and those in the yaml file 
def processYamls(
  figName,
  path=path,
  yamlNamePrefix="*",
  prefixOptions=None,# can list cellAttr etc to fine tune search 
  single=False,display=True,warningOnly=False): 
  """
  Process all jobs in a directory 
  """

  # def. 
  outCsv = path+figName+".csv"       

  # get names
  yamlNamePrefix = yamlNamePrefix.replace(".yaml","")
  if single:
      outCsv = yamlNamePrefix+"_df.csv"
      yamlNames=[yamlNamePrefix+".yaml"]
  else: 
    globTag = path+"/"+yamlNamePrefix
    if "yaml" not in globTag:
      globTag+='*yaml'

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
          columns=["trajName","tag","condVal","D","flux*A(reservoir)","flux*A(crowd)","Vol Frac","Area Frac"]
          ) 

  #print(yamlNames[0]) 
  skipped = []
  for yamlName in yamlNames:
      # open yaml
      #yamlName = path+"/"+yamlName
      # get output name 
      auxParams = GetParams(yamlName) 
      trajName = auxParams['outName']                   
      #print(trajName) 

      # compute vol frac
      volFrac = bu.CalcVolFrac(auxParams)

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
      if bu.LoadTraj(trajName,warningOnly) is None:
        skipped.append(trajName)
        continue 
      
      Di,JA, JAc, areaFrac = ProcessTraj(trajName,display=display,auxParams=auxParams) 
  
      # add to dataframe 
      df.loc[len(df.index)] = [trajName, tag, condVal,Di,JA,JAc, volFrac, areaFrac]
  
  print(df)
  print("printed %s"%outCsv)
  df.to_csv(outCsv)               

  if single:
      return True     

  if len(skipped)>0:
    print("The following were skipped ", skipped)
  

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
  msg+="  %s -fig4/-fig5/-single [yaml/dcdprefix] -all -subset [wildcard] " % (scriptName)
  msg+="""
  
 
Notes:
 It is helpful to put the subset argument in quotations 

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
      #Di,JA = processYamls("test.png",yamlName,single =True,display=True)
      #auxParams = GetParams(yamlName) 
      #ProcessTraj(yamlName,display=True,auxParams=auxParams)
      processYamls("None","./",yamlName,single=True) # need to run for printing dataframe data
      quit()

    elif(arg=="-all"): 
      processYamls("test.png",display=False,warningOnly=True,path="./")
      quit()

    elif(arg=="-subset"): 
      yamlPrefix=sys.argv[i+1]#
      processYamls("subset",display=False,warningOnly=True,path="./", yamlNamePrefix=yamlPrefix)
      quit()

  raise RuntimeError("Arguments not understood")




