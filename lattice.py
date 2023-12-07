import numpy as np
eps = 0.00   # this is used for checking sphere clashes; the effective radii are smaller than the vdw, so using zero here to allow for wiggle room
def GenerateLattice(
        nLattice,
        nRow,
        nCol=0,
        dims=[50,50],
        effDim=0.    # if positive value, will deduct this from dim so as to fit in cells (usually pass in cell diam)  
        #xThresh = None
        ):
  '''
  Generate lattice to accommodate nLattice points within the specified dims 
  '''

  dimX,dimY = dims
  if nCol==0:
    nCol = nLattice / (nRow)
  #print(nCol) 
  #nCol,nRow = 4,3
  #nLattice = nCol*nRow
  if nLattice == 1:
    nCol,nRow = 1,1

  if effDim>0:
    dimX -=effDim 
    dimY -=effDim 


  # ........######
  # .....######...
  try:
    latticeSpaceX = dimX/(nCol-1)
    latticeSpaceY = dimY/(nRow-1)
  except:
    latticeSpaceX,latticeSpaceY = dims

  #print(nCol,nRow)
  xs = np.arange(nCol)
  ys = np.arange(nRow)
  #for i in range(nCols*nRows):
  xx,yy = np.meshgrid(xs,ys)            
  xx = np.ndarray.flatten(xx)
  yy = np.ndarray.flatten(yy)
  #print(xx)
  #print(yy)
  latticePts = np.zeros([nLattice,3]) 
  latticePts[:,0] = xx*latticeSpaceX
  latticePts[:,1] = yy*latticeSpaceY

  #latticePts = np.zeros([nLattice,3]) 
  #for i in range(nLattice):
  #    try:
  #      #xi = int( np.floor(i/(nCol-1) ))
  #      xi = int( i % nCol )
  #    except:
  #      xi = 0
  #    #yi = int( i-xi*(nCol-1))
  #    yi = int( i-xi*(nCol-1))
  #    print("%2d"%i,xi,yi) 
  #    latticePts[i,0:2] = [xi*latticeSpaceX,yi*latticeSpaceY]
  

  # shift s.t. 
  if nLattice > 1: 
    latticePts[:,0] -= dimX/2.
    latticePts[:,1] -= dimY/2.
  #print(latticePts)
  

  return latticePts

def GenerateCrowderLattice(
  nCrowders,
  crowderRad=1.,
  dims=[20,20]):
  """ 
  Places crowders on a regular lattice with symmetric dimensions 
  """
  print("Placing crowders") 
  #nRow = int( np.ceil( np.sqrt( nCrowders ) )  )
  f = float(dims[1]/dims[0])  # h = f*w
  nw = np.sqrt((1/f) * nCrowders)
  nw = np.floor(nw)
  nh = np.floor(nCrowders/nw)
  nRow = int(nh)
  nCol = int(nw)
  #print(nRow,nCol)
  #print(nCrowders)
  nLattice=nRow*nCol


  #print(nRow,nLattice) 
  try: 
    width=dims[0]/(nCol-1)
    height=dims[1]/(nRow-1)
  except:
    width=1e9
  diam = 2*crowderRad

  #print(dims)
  #print(diam, width, height,nRow,nCol)
  if((width-diam)< eps  or (height-diam)< eps):
      print(diam,width,height,"d w h") 
      raise RuntimeError("Crowders are too tightly placed; check crowderRad/crowderDomain")

  #dims = [dim,dim] # square 
  #print(dims)
  latticePts = GenerateLattice(nLattice,nRow,nCol=nCol,dims=dims,effDim=diam)

  return(latticePts)


def GenerateRandomLattice( 
  crowderPosns,   # nx3 array of crowder coordinates 
  crowderRad = 10.,
  nCells = 50,
  cellRad= 1, 
  dims =  [30,30], # [um]    
  xThresh = None
  ) : 
  """ 
  Generates a randomized distribution of cells that avoids placed crowders                 
  """
  # later should adjust for asymmetric dimensions, but ignore for 
  # now 
  print("Creating cell lattice") 
  nLattice = 100
  nRow     =  int(np.sqrt(nLattice)) 
  if nCells > nLattice:
      raise RuntimeError("too many cells for lattice size") 

  latticePts = GenerateLattice(nLattice,nRow,dims=dims,effDim=(cellRad*2))

  # only keep entries to the left of xThresh
  if xThresh is not None:
    xs = latticePts[:,0]
    idx = np.where( xs <= xThresh)
    #print(np.shape(latticePts))
    latticePts = latticePts[idx[0],:]
    #print(np.shape(latticePts))

  nLattice = np.shape(latticePts)[0]
  if nCells > nLattice:
    raise RuntimeError("too many cells for lattice size %d"%nLattice)
    

  latticeIdxs= np.arange( nLattice )

  # typecast
  nCells = int(nCells)

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
  if len(allClashes) > 0:
    cellIdx = np.delete(latticeIdxs,allClashes) 
  else:
    cellIdx = latticeIdxs

  # store only those indices that are not in the violaters group 
  nCellIdx = np.shape(cellIdx)[0]
  if nCellIdx < nCells:
      raise RuntimeError("Not enough spaces to place cells") 

  
  # get random set of lattice points 

  randomized = np.random.choice( 
          np.ndarray.flatten( cellIdx ) , 
          size=nCells,
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
    crowderDims=[30,30],
    outerDims = [50,50],
    xThresh = None        
):
  #raise RuntimeError(xThresh)

  """
  Combines placement of crowders and cells onto regular lattice
  - xThresh can be used to constrain particles to the left of xThresh (so you can have assymetric starting distributions of cells) 
  """
  crowderPos = GenerateCrowderLattice(
    nCrowders, 
    crowderRad = crowderRad,
    dims=crowderDims)

  allCoords = GenerateRandomLattice( 
    crowderPosns = crowderPos, # where crowder is located 
    crowderRad = crowderRad,
    nCells = nCells,
    cellRad = cellRad,
    dims = outerDims, # [um]
    xThresh = xThresh
    ) 

  
  #print(allCoords)
  #print(crowderPos)
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




