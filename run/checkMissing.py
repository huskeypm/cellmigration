"""
Looks for jobs that haven't been processed yet
"""
import os;
import glob
files = glob.glob('*.yaml')

fileTypes=['_df.csv']
readyToProcess=[]
for fileType in fileTypes:
  for nameYAML in files:
    name=nameYAML.replace('.yaml',fileType)
    nameDCD=nameYAML.replace('.yaml',".dcd")          
    if os.path.exists(name) is not True:
      print(nameYAML+" does not have " + name )

      if os.path.exists(nameDCD):            
          readyToProcess.append(nameYAML) 

for name in readyToProcess:
    print("python ~/source//cellmigration/batchProcess.py -single %s"%name)

