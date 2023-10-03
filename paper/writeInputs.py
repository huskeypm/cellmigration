import numpy as np
'''
diam, eps, conc, +/- ATP 
crowderSize, attract


This is a script that will read a template file (cleverly
called a template.yaml) and vary a given field from 50-200%,
then print a list of yaml outputs 

You can run a bunch via listing each in -another- bash file:

python3 brown_wnonbond.py -yamlFile FILE1.yaml -run
python3 brown_wnonbond.py -yamlFile FILE2.yaml -run
....
python3 brown_wnonbond.py -yamlFile FILEN.yaml -run

Similarly, the outputs can be processed (to yield D coefficients, etc) 

python3 process.py -yamlFile FILE1.yaml -run
python3 process.py -yamlFile FILE2.yaml -run
....
python3 process.py -yamlFile FILEN.yaml -run



'''


import yaml
yamlFile='template.yaml'



def IterValues(key,nIter=3):
  dflt = auxParams[daKey] 
  # scales parameter by 0.5 to 2 
  scales = 2**np.linspace(-1,1,nIter)
  vals=dflt * scales 
  print(vals) 
  return vals 

def WriteYaml(contents, fileName):
  # write fitness and best params to file
  with open (fileName, 'w') as file:
    documents = yaml.dump(contents, file)
  print("Created yaml file ",fileName)

##
##
##

with open(yamlFile, 'r') as file:
  auxParams = yaml.safe_load(file)

for key in auxParams.keys():
  print(key,auxParams[key])

keys=['nParticles']

daKey = keys[0]

vals = IterValues(daKey)

for val in vals:
  val = float(val) 
  #print(val) 
  outParams = auxParams.copy()
  outParams[daKey] = val
  keyName="_%s%f"%(daKey,val)
  outParams['outName']=auxParams['outName']+keyName
  writeName = "out"+keyName+".yaml"


  yaml.safe_dump(outParams, sort_keys=False)
  WriteYaml(outParams,writeName)
  
