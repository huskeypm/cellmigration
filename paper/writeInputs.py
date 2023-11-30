import numpy as np
'''


This is a script that will read a template file (cleverly
called a template.yaml) and vary a given field from 50-200%,
then print a list of yaml outputs 

* You can run a bunch via listing each in -another- bash file:

python3 brown_wnonbond.py -yamlFile FILE1.yaml -run
python3 brown_wnonbond.py -yamlFile FILE2.yaml -run
....
python3 brown_wnonbond.py -yamlFile FILEN.yaml -run

* The jobs can be processed en masse using 
batchProcess.py 

* Or, the outputs can be processed (to yield D coefficients, etc) 

python3 process.py -yamlFile FILE1.yaml -run
python3 process.py -yamlFile FILE2.yaml -run
....
python3 process.py -yamlFile FILEN.yaml -run

* To implement
- diam [cellRad | crowderRad?] 
- eps [yScale/yDiam - tricky, since attraction needs to be applied i mjust a narrow region; might be easier w crowders] 
- conc [nParticles, DONE] 
- +/- ATP [xScale] 
- attract [cellAttr, crowderAttr]

'''


##
## PARAMS 
##
import yaml

varIter=5 # range of key values
varRuns=3    # number of time each condition is run 

class Params():
    def __init__(self,
      cmd = "python3 ../brown_wnonbond.py -yamlFile ",
      path="/home/pkh-lab-shared/migration/test/"   # path for outfiles 
            ):

      self.cmd = cmd
      self.path=path

params=Params()

##
## Func 
##

def RescaleValues(dflt,nIter=3,key="key"):
  # scales parameter by 0.5 to 2 
  scales = 2**np.linspace(-1,1,nIter)
  vals=dflt * scales 
  print(key,vals) 
  return vals 

def WriteYaml(contents, fileName,verbose=True):
  # write fitness and best params to file
  with open (fileName, 'w') as file:
    documents = yaml.dump(contents, file)
  if verbose:
    print("Created yaml file ",fileName)

# for a given key, write rescaled calues 
def WriteIterativeYaml(auxParams,daKey,varIter=3,varRuns=2,writeFile=None):
  #print('path',params.path,'exec',params.cmd) 

  dflt = auxParams[daKey] 
  vals = RescaleValues(dflt,nIter=varIter,key=daKey) 

  # range of vals
  cmds = []
  for val in vals:
    val = float(val) 
    #print(val) 
    outParams = auxParams.copy()
  
    # give a tag to determine which value was iterated
    outParams['tag']=daKey
    outParams[daKey] = val
   
    # replicates per val 
    for run in range(varRuns):
      keyName="_%s%f_%.2d"%(daKey,val,run)
      #outParams['outName']=path+"/"+auxParams['outName']+keyName
      #writeName = outParams['outName']+".yaml"
      outParams['outName']=auxParams['outName']+keyName
      writeName = outParams['outName']+".yaml"
      fullWriteName = params.path+"/"+writeName

      yaml.safe_dump(outParams, sort_keys=False)
      WriteYaml(outParams,fullWriteName,verbose=False)
      daCmd = (params.cmd+" -yamlFile " + writeName+" -run")
      cmds.append(daCmd)
  ##########################     
  if writeFile is not None:
    with open(writeFile,'w') as f:
      [f.write(i+"\n") for i in cmds]
      print("Wrote "+writeFile)
  else:
      [print(i) for i in cmds]

    
#!/usr/bin/env python
##
## MAIN 
##
def main(yamlFile, keys=['nParticles']):

  varIter = 3
  with open(yamlFile, 'r') as file:
    auxParams = yaml.safe_load(file)
  
  for key in auxParams.keys():
    print(key,auxParams[key])
  
  daKey = keys[0]
  for daKey in keys:
    WriteIterativeYaml(auxParams,daKey,varIter=varIter) 
  
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

  yamlFileNoCrwd='template.yaml'
  yamlFileNoCrwdATP='template_atp.yaml'
  yamlFileWCrwd='template_crowders.yaml'
  yamlFileWCrwdATP='template_crowders_atp.yaml'
  keys=['nParticles']

  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-fig4"):
      main(yamlFile=yamlFileNoCrwd, keys=['nParticles','cellRad','cellAttr'])
      main(yamlFile=yamlFileNoCrwdATP, keys=['nParticles','cellRad','cellAttr'])
      quit()
    elif(arg=="-fig5"):
      #main(yamlFile=yamlFileWCrwd, keys=['nCrowders','crowderRad','xScale'])
      main(yamlFile=yamlFileWCrwd, keys=['nCrowders','crowderRad','crowderAttr'])
      main(yamlFile=yamlFileWCrwdATP, keys=['nCrowders','crowderRad','crowderAttr'])
      quit()
    elif(arg=="-fig6"):
      main(yamlFile=yamlFileWCrwd, keys=['yScale','ySize'])                      
      print("yscale stuff needs implemented; mainly to vary width of channel like dumbbell shape") 
      print("if yScale|ySize, yPotent=True") 
      print("consder just using crowders" ) 
      quit()
  





  raise RuntimeError("Arguments not understood")




