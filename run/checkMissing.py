import os;
import glob
files = glob.glob('*.yaml')

fileTypes=['.dcd','_df.csv']
fileTypes=['_df.csv']
for fileType in fileTypes:

  for file in files:
    name=file.replace('.yaml',fileType)
    if os.path.exists(name) is not True:
      print(file+" does not have " + name )


