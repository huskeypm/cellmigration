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



def ProcessTraj(caseName,display=False): 
    # LoadTraj
    caseName = caseName.replace('.yaml',"")
    traj = bu.LoadTraj(caseName)

    #, get RDF
    bu.CalcRDFs(traj,"RC","AC")
    #mask1='@RC'
    #mask2=':21@AC'# last atom (crowder) 
    #bu.CalcRDF(traj,mask1=mask1,mask2=mask2)  
#
    ## get J,D    
    Di=bu.CalcD(traj,mask='@RC',csvName=caseName)                        
    xThresh = 300 # AA
    xThresh = 600 # AA
    JA=bu.CalcFlux(traj,mask='@RC',display=display,xThresh=xThresh,caseName=caseName)
    print("Di %f J %f"%(Di,JA))

    # get 2D histogram of populations
    bu.CalcProbDist(traj,mask='@RC',caseName=caseName,display=display)

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
  yamlNamePrefix = yamlNamePrefix.replace(".yaml","")
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
  skipped = []
  for yamlName in yamlNames:
      # open yaml
      #yamlName = path+"/"+yamlName
      with open(yamlName, 'r') as file:                                
        auxParams = yaml.safe_load(file)
      # get output name 
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
      
      Di,JA = ProcessTraj(trajName,display=display) 
  
      # add to dataframe 
      df.loc[len(df.index)] = [trajName, tag, condVal,Di,JA,volFrac]
  
  if single:
      return Di,JA
  
  print(df)
  outCsv = path+figName+".csv"       
  print("printed %s"%outCsv)
  df.to_csv(outCsv)               

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
      #Di,JA = processYamls("test.png",yamlName,single =True,display=True)
      Di,JA = ProcessTraj(yamlName,display=True)
      quit()

    elif(arg=="-all"): 
      processYamls("test.png",display=False,warningOnly=True,path="./")
      quit()

  raise RuntimeError("Arguments not understood")




