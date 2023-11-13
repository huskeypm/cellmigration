import numpy as np



def GenerateLattice(
        nLattice,
        nRow,
        dim,
        effDim=0.    # if positive value, will deduct this from dim so as to fit in particles (usually pass in particle diam)  
        ):

  if effDim>0:
      dim-=effDim 

  try:
    latticeSpace = dim/(nRow-1)
  except:
    latticeSpace = dim


  
  latticePts = np.zeros([nLattice,3]) 
  for i in range(nLattice):
      xi = int( np.floor(i/nRow) )
      yi = int( i-xi*nRow )
      latticePts[i,0:2] = [xi*latticeSpace,yi*latticeSpace]
  

  # shift s.t. 
  latticePts[:,0] -= dim/2.
  latticePts[:,1] -= dim/2.

  return latticePts

def GenerateCrowderLattice(
  nParticles,
  crowderRad=1.,
  dim=20):
  """ 
  Places crowders on a regular lattice
  """
  nRow = int( np.ceil( np.sqrt( nParticles ) )  )
  nLattice=nRow**2
  try: 
    width=dim/(nRow-1)
  except:
    width=1e9
  diam = 2*crowderRad

  if(diam>(width-0.01)):
      raise RuntimeError("Crowders are too tightly placed; check crowderRad/crowderDomain")

  latticePts = GenerateLattice(nLattice,nRow,dim,effDim=diam)

  return(latticePts)


def GenerateRandomLattice( 
  # PKH turn into nCrowder x 3 array 
  #crowderPos = np.array([0,0,0]), # where crowder is located 
  crowderPosns,   # nx3 array of crowder coordinates 
  crowderRad = 10.,
  nParticles = 50,
  particleRad= 1, 
  dim =  30 # [um]    
  ) : 
  """ 
  Generates a distribution of cells that avoids placed crowders                 
  """
  nLattice = 100
  nRow     =  10 
  latticePts = GenerateLattice(nLattice,nRow,dim,effDim=(particleRad*2))
  latticeIdxs= np.arange( np.shape(latticePts)[0])

  # typecast
  nParticles = int(nParticles)

  # find crowder positions that conflict w lattice 
  # PKH iterate over each crowder position 
  nCrowders = np.shape(crowderPosns)[0]
  allClashes = []               
  for i in range(nCrowders):
    minDistSqd = (latticePts[:,0] - crowderPosns[i,0])**2
    minDistSqd+= (latticePts[:,1] - crowderPosns[i,1])**2
    #crowderIdx = np.argmin(minDist)
    #print(crowderIdx)
  
    # 'remove' conflicting points 
    # PKH find all cells that VIOLATE crowderRad
    #cellIdx = np.argwhere(minDist > crowderRad**2) 
    clashIdx = np.argwhere(minDistSqd <= crowderRad**2) 
    #clashIdx = np.ndarray.flatten(clashIdx)
    #print(clashIdx)
    if len(clashIdx)>0:
      for clash in clashIdx: 
        allClashes.append(clash)

  # remove dupes 
  allClashes = np.asarray(allClashes)
  allClashes = np.ndarray.flatten( allClashes) 
  allClashes = np.unique(allClashes)

  # remove conflicting entries 
  cellIdx = np.delete(latticeIdxs,allClashes) 
  latticeIdxs = cellIdx

  # store only those indices that are not in the violaters group 

  nCellIdx = np.shape(cellIdx)[0]
  if nCellIdx < nParticles:
      raise RuntimeError("Not enough spaces to place cells") 

  
  # get random set of lattice points 
  #print(np.shape(cellIdx)) 

  randomized = np.random.choice( 
          np.ndarray.flatten( cellIdx ) , 
          size=nParticles,
          replace=False) # dont reuse lattice points 

  #for i in cellIdx:
  #    print(latticePts[i,0:2], minDist[i])
  #print( latticePts[cellIdx,:] ) 
  #print( minDist[cellIdx])
  #print( np.shape(cellIdx) ) 
  
  cellCoords = latticePts[randomized,:]
  return cellCoords

def GenerateCrowdedLattice(
    nCrowders, 
    nCells,
    crowderRad=10.,
    cellRad=1.,
    crowdedDim=30,
    outerDim = 50):

  """
  Combines placement of crowders and cells onto regular lattice
  """
  crowderPos = GenerateCrowderLattice(
    nCrowders, 
    crowderRad = crowderRad,
    dim=crowdedDim)

  #print(outerDim) 
  allCoords = GenerateRandomLattice( 
    crowderPosns = crowderPos, # where crowder is located 
    crowderRad = crowderRad,
    nParticles = nCells,
    particleRad = cellRad,
    dim = outerDim # [um]
    ) 

  
  return crowderPos, allCoords 

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

  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-generate"):
      GenerateCrowdedLattice(16,20)
      quit()

  





  raise RuntimeError("Arguments not understood")




