#!/usr/bin/env python
import numpy as np
import brown_util as bu
import yaml 


def processPKL(params):
  #print(params.keys())
  #trajOutName = "./tests/crowder.pkl"
  trajOutName=params['outName']+".pkl"
  
  ts,xs,ys, nUpdates, nParticles = bu.LoadPKLData(trajOutName)
  msds = bu.meanSquareDisplacements(xs,ys,nUpdates)
  texp, msdfit,D=bu.CalcMSD(ts,msds)
  D=D[0]
  return D 

#!/usr/bin/env python
import sys
##################################
#
# Revisions
#       10.08.10 inception
#
##################################

#
# ROUTINE  
#
def doit(yamlFile):
  # read
  with open(yamlFile, 'r') as file:
    auxParams = yaml.safe_load(file)
  file.close()

  # get D 
  D = processPKL(auxParams)

  # write 
  outParams = auxParams.copy()
  outParams['D'] = "%f"%D
  outYamlFile = yamlFile.replace(".yaml","_out.yaml")
  with open(outYamlFile, 'w') as file:
    yaml.safe_dump(outParams,file)                   
  file.close()

  # display 
  daStr=""
  for key in outParams.keys():
      daStr+=key+","
      daStr+="%s"%outParams[key]+","
  print(daStr)    


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

  # Loops ovyer each argument in the command line 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-yaml"):
      arg1=sys.argv[i+1] 
      doit(arg1)
      quit()
  





  raise RuntimeError("Arguments not understood")




