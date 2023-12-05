
import numpy as np 
#from numpy import random
import random 
from scipy.spatial.distance import cdist, pdist, squareform 
def UpdateBoundary(simulation,paramDict,x,nCells):
  """
  Apply absorbing bounadary on right hand side. Deletes particles on right; move them to left
  x - current atom positions
  nCells
  """
  # find things near other boundary 
  xys = x[:nCells,0:2]
  absorbing = paramDict["domainXDim"]/2 - paramDict["absorbingMargin"]

  #absorbing = -26  
  #print(xys, absorbing)
  moved = np.where( xys[:,0] > absorbing )
  nMoved = len(moved[0])
  if nMoved < 1: 
    return 0 
  print("moving %d particle (%f)"%(nMoved,absorbing)) 

  # define interval of lattice points that can accommodate moved particles 
  cellDiam= 2*paramDict["cellRad"] * 0.9 # cells usually fit kind of tight
  #s1 = np.array([(0,0), (0,1), (1,0), (1,1)])
  daxMin = -paramDict["domainXDim"]/2 + 0.5*cellDiam 
  daxMax = daxMin+paramDict["absorbingMargin"]
  dayMin = -paramDict["domainYDim"]/2 + 0.5*cellDiam 
  dayMax = -dayMin
  xvalues = np.linspace( daxMin,daxMax, int((daxMax-daxMin)/cellDiam) + 1) 
  yvalues = np.linspace( dayMin,dayMax, int((dayMax-dayMin)/cellDiam) + 1) 
  xx, yy = np.meshgrid(xvalues, yvalues)
  trial = np.dstack([xx, yy]).reshape(-1, 2)
  #print(trial) 

  # check distance between lattice points and existing cells 
  #print(cdist(xys,trial))
  v=cdist(xys,trial).min(axis=0)
  #print(v) # distance to closest cell for each lattice pint
  #print(distance.cdist(s1,s2).argmin(axis=0))
  args = np.where(v > cellDiam)[0]
  try:
    z = random.sample(set(args), nMoved)
  except:
    print(xys)
    print(trial) 
    print(v) 
    raise RuntimeError("not enough room to accommodate particle (%d)"%nMoved) 
  #print('z',z)
  #print(trial[z])

  # place edge particles on selected spots 
  xyz[z,:] = trial[z,:] 
  x[:,0:2]= xyz[:,]
  simulation.context.setPositions(x)

  return nMoved 
  
  

def UpdateStates(csi,x,t,idxsCells,idxsA,idxsB,paramDict):      
      """
      find close contacts between cells and crowders 
      this seems inefficient; seems like openmm should have something 
      should package this into brown_util 

      csi - CellSystem object
      x - current coordinates
      t - current time 
      idxsCells - indices of cells
      """
      dists=pdist(x)
      dists = squareform(dists)
      daMax = np.max(dists)
      #for i in range(10): # easy way to prevent self interction
      #  dists[i,i]=daMax
      # for debug
      #dists = np.array(dists,int)
      #print(dists)
      closeA=bu.GetContacts(dists,idxsA,thresh=paramDict["contactDistCrowderA"])
      closeB=bu.GetContacts(dists,idxsB,thresh=paramDict["contactDistCrowderB"])

      # update contacts list in state object
      #csi.UpdateContactsA( np.ones(nCells) )
      #csi.UpdateContactsB( np.ones(nCells) )
      csi.UpdateContactsA( closeA[idxsCells] )              
      csi.UpdateContactsB( closeB[idxsCells] )              

      # determine next state based
      #csi.PrintStates()
      csi.EvalTransitions(t)#,nCells,cs.stateMatrix)
      #csi.PrintStates()
      stateSums = csi.SumStates()

      closeAs = np.sum(closeA[idxsCells])
      closeBs = np.sum(closeB[idxsCells])

      return stateSums,closeAs,closeBs


