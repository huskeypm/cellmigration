#!/cm/shared/apps/anaconda3/bin/python3
# adapted from https://gist.github.com/rmcgibbo/6094172
# code also adapted from runner.py written by Ben Chun

# v4 includes support for states via gillespie  
#
import matplotlib.pylab as plt
from scipy.spatial.distance import cdist,pdist, squareform
import random


"""
Propagating 2D dynamics on arbitrary potential using OpenMM.
Currently, we just put a harmonic restraint on the z coordinate,
since OpenMM needs to work in 3D. This isn't really a big deal, except
that it affects the meaning of the temperature and kinetic energy. So
take the meaning of those numbers with a grain of salt.
"""

import socket
isKant = False
if socket.gethostname()=='kant':
    isKant=True

# if kant installastion
if isKant:
  from simtk.unit import kelvin, picosecond, femtosecond, nanometer, dalton
  import simtk.openmm as mm
  from simtk.openmm.app import *                    

# otherwise 
else: 
  import openmm as mm
  from openmm.unit import kelvin, picosecond, femtosecond, nanometer, dalton
  from openmm.app import * # PDBFile, DCDReporter

import lattice 
import brown_util as bu
import states
import numpy as np

## INIT 
import platform as pf
if pf.system()=='Darwin' or isKant:
  print("Running on mac; assuming no CUDA")
  platform = mm.Platform.getPlatformByName('CPU')
  properties = {}
else:
  print("Linux") 
  platform = mm.Platform.getPlatformByName('CUDA')
  properties = {'Precision': 'double'}


min_per_hour = 60  #


# https://demonstrations.wolfram.com/TrajectoriesOnTheMullerBrownPotentialEnergySurface/#more
# adding external potential to keep things within a cylinder, based on this example
# http://docs.openmm.org/latest/userguide/application/04_advanced_sim_examples.html 
class CustomForce(mm.CustomExternalForce):
    """OpenMM custom force for propagation on the Muller Potential. Also
    includes pure python evaluation of the potential energy surface so that
    you can do some plotting"""
    aa = [5e-3]       
    bba= [1]   # power for x (1 - linear) 
    bb = [5e-3]
    YY = [0]
    # z - potential is defined below 

    # for square potential 
    MINX=[-50e3]  # updated below 
    MAXX=[ 50e3]
    MINY=[-10e3]
    MAXY=[ 10e3]



    def __init__(self,
        paramDict 
        ):
        
        pD = paramDict
        yPotential = pD["yPotential"]
        xPotential = pD["xPotential"]
        containmentPotential = pD["containmentPotential"]
        self.MINX[0]=-pD["domainXDim"]/2.
        self.MAXX[0]= pD["domainXDim"]/2.
        self.MINY[0]=-pD["domainYDim"]/2. # 4.
        self.MAXY[0]= pD["domainYDim"]/2. # 4.
        self.YY[0] = -pD["domainYDim"]/2.

        # start with a harmonic restraint on the Z coordinate
        expression = '100.0 * z^2'

 
        # chemoattractant gradient ; assume RHS is maximal ATP
        # s = ( c(x=max)-c(x=min) ) /(max-min) ==> s = c/(2*xmax) 
        # b = 1/2 c
        # y = x * s + b ==> y = x c/(2m) + 1/2 (m*c/m) 
        # ==> y = c/2m * (x + m)
        ac    = pD["chemoAttr"]/(2*self.MAXX[0])
        #print(ac)
        self.XX = [-self.MAXX[0]] 
        #print("TODO get rid of indexing a[0]") 
        self.aa[0] = -1 * pD["xScale"] * ac # make attractive for U = -xScale * c       


        # y 
        self.bb[0] = -1 *  pD["yScale"]
         

        # any changes here must be made in potential() below too 
        j=0   # TODO superfluous, remove later 

        if yPotential is False:
          self.bb[j]=0.
        if xPotential is False:
          self.aa[j]=0.
             
        # add the terms for the X and Y
        fmt = dict(
                aa=self.aa[j], XX=self.XX[j], bb=self.bb[j], bba=self.bba[j], YY=self.YY[j],
                MINX=self.MINX[j], MAXX=self.MAXX[j],MINY=self.MINY[j], MAXY=self.MAXY[j], # for square
                )

        # y parabola 
        #expression += '''+ {bb} * (y - {YY})^4'''.format(**fmt)
        expression += '''+ {bb} * exp(-(y - {YY})/25.)'''.format(**fmt)

        # xExpression / gradient 
        #expression += '''+ {aa} * (x - {XX})^4'''.format(**fmt)
        # my potential for the x direction
        # xPotential = aa*(xParticle-x0)      <--- if bba=1
        expression += '''+ {aa} * (x - {XX})^{bba}  '''.format(**fmt)
        expression +=" "


        # I bet this can be converted using sympy or something, but for now just use this 
        # plotme 
        def func(x,y):
            z=0
            #v =100.0 * z**2+ 0.0 * (y - 0)**4+ 1.0 * (x - -50.0)**1   # +10*(max(0, -50.0-x) + max(0, x-50.0) + max(0, -25.0-y) + max(0, y-25.0));
            v=100.0 * z**2  + -1.0 * np.exp(-(y - -25.0)/25.)+ -0.01 * (x - -50.0)**1
            return v
        vfunc = np.vectorize(func)
        display = False 
        if display:
          #print(vfunc(1,2))
          #xx,yy = np.meshgrid(0:10,0:10)
          xx,yy = np.mgrid[-pD['domainXDim']/2:pD['domainXDim']/2:1,-pD['domainYDim']/2:pD['domainYDim']/2:1]
          V = vfunc(xx,yy)

          plt.pcolormesh(xx,yy,V,shading='gouraud')
          plt.axis('equal')
          plt.colorbar()
          plt.gcf().savefig("energy.png")
          quit()


        # cylindrical container in xy plane 
        if containmentPotential=='cylinder':
          print('adding cylindrical containment potential') 
          expression+='+10*max(0, r-10)^2; r=sqrt(x*x+y*y+z*z)'

        elif containmentPotential=='square':
          print('adding square containment potential') 
          # works 
          #expression+='+10*(max(0, -10-x) + max(0, x-10));'
          expression+='''+10*(max(0, {MINX}-x) + max(0, x-{MAXX}) + max(0, {MINY}-y) + max(0, y-{MAXY}));'''.format(**fmt)
          #expression+="+100*(0*step(10-x) + (1-step(-10-x)))"#+
#                             step(10-y) + step(-10-y))     "

        print(expression)
                               
        super(CustomForce, self).__init__(expression)


# allocate instance 
import parameters 
params = parameters.Params()

print("MOVE UPDATES TO BROWN_UTIL") 
from scipy.spatial import distance
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
  cellDiam= 2*paramDict["cellRad"] * 1.0 # to allow some wiggle room 
  #s1 = np.array([(0,0), (0,1), (1,0), (1,1)])
  daMax = paramDict["domainYDim"]/2 - 0.5*cellDiam 
  #print(daMax)
  ys = np.linspace(-daMax,daMax, int(2*daMax/cellDiam + 1)) 
  trial = np.zeros([np.shape(ys)[0],2])
  trial[:,0] = -paramDict["domainXDim"]/2 + 0.5*cellDiam
  trial[:,1] = ys
  #print(trial) 

  # check distance between lattice points and existing cells 
  #print(cdist(xys,trial))
  v=cdist(xys,trial).min(axis=0)
  #print(v) # distance to closest cell for each lattice pint
  #print(distance.cdist(s1,s2).argmin(axis=0))
  args = np.where(v > cellDiam)[0]
  #print(args)
  try:
    z = random.sample(set(args), nMoved)
  except:
    raise RuntimeError("not enough room to accommodate particle") 
  #print('z',z)
  #print(trial[z])

  # place edge particles on left 
  #print('to move',moved) 
  #print(xys)
  xys[moved] = trial[z]                     
  x[:nCells,0:2]=xys
  #print(xys)
  check = squareform(pdist(xys))
  #print(check)
  #print(check[moved[0]])

  # update positions 
  simulation.context.setPositions(x)                    
  #print(i,np.min(xs),np.max(xs))

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




import yaml
def runBD(
  # each particle is totally independent, propagating under the same potential
  display=False,
  yamlFile=None
  ): 
  paramDict = params.paramDict

  if yamlFile is None:
    print("YAML file was not provided - using default parameters.") 
   
  else:
    # load yaml
    with open(yamlFile, 'r') as file:
      auxParams = yaml.safe_load(file)
    
    for key in auxParams.keys():
      # check if key is defined already (usually bad if not) 
      if key not in paramDict.keys():
          raise RuntimeError(key+" is not defined") 

      # assign 
      paramDict[key] = auxParams[key]

      if "trajOutName" in key:
          raise RuntimeError(key+" is now deprecated. Use outName instead")
      print("Adding %s="%(key) , auxParams[key])

  params.update()
    
  # place particles 
  # TODO: start w preequilibrated box or get cells from expt 
  nCells = int(paramDict["nCells"])
  nCrowders = int(paramDict["nCrowders"])


  if paramDict["effectiveRad"] is None:
    crowderRad = paramDict['crowderRad'] 
  else:
    crowderRad = paramDict['effectiveRad'] 

  outerDims=[paramDict["domainXDim"], paramDict["domainYDim"]]
  crowderPos, cellPos = lattice.GenerateCrowdedLattice(
          nCrowders,nCells,
          crowderRad,paramDict['cellRad'],
          crowdedDim=paramDict["crowderDim"], # [um] dimensions of domain containing crowders (square)  
          outerDims=outerDims,
          xThresh=paramDict["xThresh"]
          )  # generate crowders

  newCrowderPos = np.shape(crowderPos)[0]
  if (newCrowderPos != nCrowders):
    print("WARNING: increasing nCrowders to the nearest 'square-rootable' value: %d"
            %newCrowderPos)
    nCrowders = newCrowderPos 

  # cells first, then crowders
  nTot = nCells + nCrowders
  startingPositions = np.concatenate((cellPos,crowderPos))
  #print(np.shape(crowderPos))
  #print(np.shape(cellPos))
  #print(np.shape(startingPositions))

  # align everything to xy plane 
  startingPositions[:,2] = 0.


  print("WARNING: adding random noise, since there's a weird bug with particles placed at origin 0,0,0")
  print("TODO try removing") 
  startingPositions[:,0]+=1e-3*np.random.rand(nTot)
  startingPositions[:,1]+=1e-3*np.random.rand(nTot)
  ###############################################################################

  
  system = mm.System()

  ## define outputs for coordinates
  #import calculator as calc 
  trajOutPfx=paramDict["outName"]
  trajOutName = trajOutPfx+".pkl"
  pdbFileName = trajOutPfx+".pdb"
  dcdFileName = trajOutPfx+".dcd"
  # define arbitrary pdb
  nm_to_Ang=10
  sp_Ang = startingPositions*nm_to_Ang # default is nm in program, but pdb/dcd use Ang     
  bu.genPDBWrapper(pdbFileName,nCells,nCrowders,sp_Ang)
  #calc.genPDBWrapper(pdbFileName,nTot,startingPositions)
  # add to openmm
  pdb = PDBFile(pdbFileName) 

  # Configure dcd                    
  dumpSize = paramDict['dumpSize'] # 100 
  dcdReporter = DCDReporter(dcdFileName, dumpSize)

  print("ADD ME")
  #simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
  #      potentialEnergy=True, temperature=True))


  # define external force acting on particle 
  # WARNING: ALWAYS treat particles (cells) before crowders, since that's how they're printed in the pdb file 
  # better fix for this????
  customforce = CustomForce(paramDict)
  cfi=0
  for i in range(nCells):      
      system.addParticle(paramDict["mass"])
      customforce.addParticle(cfi, [])
      cfi+=1
  for i in range(nCrowders):      
      #system.addParticle(paramDict["mass"]*1e4)
      system.addParticle(0.)    # enforces them to be fixed 
      customforce.addParticle(cfi, [])
      cfi+=1
  system.addForce(customforce) # <-- PKH should this be added earlier to keep things in zi (i think it can be) 


  
  nonbond = mm.NonbondedForce()
  for i in range(nCells):      
     # chg, LJ sigma (nm), well-depth (kJ)
    [q,d,w] = paramDict["cellChg"], paramDict['cellRad'], paramDict['cellAttr']
    nonbond.addParticle(q,d,w)                
  for i in range(nCrowders):      
    [q,d,w] = paramDict["crowderChg"], paramDict['crowderRad'], paramDict['crowderAttr']
    nonbond.addParticle(q,d,w)                
  # Add force to the System
  system.addForce(nonbond)

  #integrator = mm.LangevinIntegrator(temperature, friction, timestep)
  integrator = mm.BrownianIntegrator(paramDict["temperature"], paramDict["friction"], paramDict["timestep"])
  simulation = Simulation(pdb.topology, system,integrator,platform,properties) 
  
  simulation.context.setPositions(startingPositions)
  simulation.context.setVelocitiesToTemperature(paramDict["temperature"])

  # dcd writing
  simulation.reporters.append(dcdReporter)
  
  nUpdates = int(paramDict["nUpdates"])
  totTime = nUpdates *  paramDict["frameRate"]  # [1 min/update]

  ts = np.arange(nUpdates)/float(nUpdates) * totTime             
  xs = np.reshape( np.zeros( nTot*nUpdates ), [nTot,nUpdates])
  ys = np.zeros_like(xs)
  

  # minimize to reconcile bad contacts
  minimize = True   
  if minimize:
    print("Minimizing system") 
    print("If stalls, there's probably a big clash somewhere") 
    simulation.minimizeEnergy() # don't do this, since it will move the crowders too 
  else:
    print("WARNING: not minimizing; should have per-crowder constraints") 

  iters = nUpdates/10.
  import timeit
  start = timeit.default_timer()
  diffMax=0.


  #state.init(nCell)
  #idxCells    = 0  # range 
  #idxCrowderA = 0+nCell
  #idxCrowderB = indexCrowderA +nCrowder/2
  nB = int(nCrowders/2)
  nA = nCrowders-nB
  idxsCells = np.arange(nCells) 
  idxsA = np.arange(nA)+nCells
  idxsB = np.arange(nB)+nCells+nA
  #print(idxsCells,idxsA,idxsB)

  #
  # START ITERATOR 
  #
  print("Running dynamics") 
  csi = states.CellSystem(nCells=nCells)
  csi.UpdateK01(paramDict["K01"])
  csi.UpdateK12(paramDict["K12"])
  csi.UpdateK20(paramDict["K20"])
  updateStates = 10
  stateUpdates = int(nUpdates/updateStates)
  closeAs = np.zeros(stateUpdates)
  closeBs = np.zeros(stateUpdates)
  stateSums= []
  movedParticles = 0 

  for i in range(nUpdates): 
      if (i % iters)==0:
        print("(%d/10)..."%(i/iters+1))

      # get positions at ith cycle 
      x = simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(nanometer)

      # check for large displacement 
      try:
        diff =np.abs(x-xprev)
        diff= np.max(diff)
        #print(diff)
      except:
        diff=0.
      diffMax = np.max([diffMax,diff])
      if diffMax > 100:
          print("Something happened at step %i (large displacement); stopping"%i)
          break    

      # should only use when asymmetric distro is used 
      if paramDict["absorbingBoundary"]:
        mvd = UpdateBoundary(simulation,paramDict,x,nCells)
        movedParticles+=mvd
        # probably want to double check that velocities aren't reset 
        


      if (paramDict['states'] and (i % updateStates) == 0):
        t=i  # time/not sure what units/values to use yet ??????
        j = int(i/updateStates)
        stateSum,closeA, closeB = UpdateStates(csi,x,t,idxsCells,idxsA,idxsB,paramDict)
        closeAs[j] = closeA
        closeBs[j] = closeB
        stateSums.append(stateSum) 


      # get particle's positions 
      #xsi[:,i] = x[:,0]*nm_to_Ang
      #ysi[:,i] = x[:,1]*nm_to_Ang
      #j=1 # particle 2
      #print(xsi[j,i],ysi[j,i])
      #print("coord",xs[10:,i],ys[10:,i])
      #print("------") 
      #print(x[0:5,:])
      
      # integrate 
      #integrator.step(paramDict["nInteg"]) # 100 fs <-->  
      simulation.step( paramDict["nInteg"] ) 



      # 
      xprev = np.copy(x) 
  #
  # END ITERATOR 
  #
  print("Moved %d particles\n"%movedParticles)
  if diffMax > np.min([paramDict["crowderRad"],paramDict["cellRad"]]):
      print("Diff max %f"%diffMax)
      print("WARNING: step size is greater than particle sizes, which may give bad results;increase friction") 
  stop = timeit.default_timer()
  print('Time(s): ', int(stop - start))

  # plot states 
  if paramDict['states']:
    plt.figure()
    plt.title("state simulation") 
    fig, axl = plt.subplots()
    ts = np.arange(np.shape(closeAs)[0])
    axl.scatter(ts,closeAs,label='#cells contacting A')
    axl.scatter(ts,closeBs,label='B')
    axl.legend(loc=1)
    axr = axl.twinx()
    daKeys = stateSums[0].keys()
    for key in daKeys:
        l = [i[key] for i in stateSums]
        axr.plot(l,label=key)
        axr.legend(loc=2)
    plt.gcf().savefig("close.png") 

  #if display:
  #    plt.show() 

  # package data 
  printPKL=False 
  if printPKL:
    ar = [ts,xs,ys]
    import pickle as pkl
    if trajOutName is not None:
      if "pkl" not in trajOutName:
        trajOutName+=".pkl"
      file = open(trajOutName, 'wb') 
      pkl.dump(ar,file)        
      file.close()
    print("WARNING: pdb, dcd, and pkl are reported in Angstrom, while states are in nanometer!!")
  
  return ts,xs, ys 

     

#!/usr/bin/env python
import sys

# prints variables availasble for yaml
def printVariables():
    print("These variables can be added to the yaml file")
    params = Params()

    paramDict = params.paramDict
    for key, value in paramDict.items():
        print("key/val ",key,value) 




#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -validation/-printVar" % (scriptName)
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


  if len(sys.argv) < 2:
      raise RuntimeError(msg)



  #fileIn= sys.argv[1]
  #if(len(sys.argv)==3):
  #  1
  #  #print "arg"

  display=False 
  yamlFile = None 
  
  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-display"):
      display=True
    if(arg=="-yamlFile"):
      yamlFile= sys.argv[i+1]
    if(arg=="-outName"):
      params.paramDict["outName"]=sys.argv[i+1]

    #
    # Run modes 
    # 
    if(arg=="-validation"):
      #arg1=sys.argv[i+1] 
      runBD(display=display,yamlFile=yamlFile)
      quit()

    if(arg=="-run"):
      runBD(display=display,yamlFile=yamlFile)
      quit()

    if(arg=="-printVar"):
      printVariables()
      quit()
  





  raise RuntimeError("Arguments not understood")




