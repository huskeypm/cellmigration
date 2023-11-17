#!/cm/shared/apps/anaconda3/bin/python3
# adapted from https://gist.github.com/rmcgibbo/6094172
# code also adapted from runner.py written by Ben Chun

# Not convinced that the LJ terms are working appropriately, therefore 
# I'm troublshooting based on the helloargon openmm example 
# http://docs.openmm.org/latest/userguide/library/03_tutorials.html#helloargon-program

"""Propagating 2D dynamics on the muller potential using OpenMM.
Currently, we just put a harmonic restraint on the z coordinate,
since OpenMM needs to work in 3D. This isn't really a big deal, except
that it affects the meaning of the temperature and kinetic energy. So
take the meaning of those numbers with a grain of salt.
"""

#from simtk.unit import kelvin, picosecond, femtosecond, nanometer, dalton
# if kant installastion
#import simtk.openmm as mm
#from simtk.openmm.app import *                    
# otherwise 
import openmm as mm
from openmm.unit import kelvin, picosecond, femtosecond, nanometer, dalton
from openmm.app import * # PDBFile, DCDReporter
import lattice 
import brown_util as bu
import numpy as np

## INIT 
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
    XX = [0] # need to adjust this to the domain size 
    YY = [0]
    # z - potential is defined below 

    # for square potential 
    MINX=[-50e3]
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
        self.MINX[0]=-pD["domainDim"]/2.
        self.MAXX[0]= pD["domainDim"]/2.
        self.MINY[0]=-pD["domainDim"]/2. # 4.
        self.MAXY[0]= pD["domainDim"]/2. # 4.

        # start with a harmonic restraint on the Z coordinate
        expression = '100.0 * z^2'

 
        # chemoattractant gradient ; assume RHS is maximal ATP
        # c(x=len) = chemoAttr * ac *x ---> ac = c/lenDom 
        ac    = pD["chemoAttr"]/self.MAXX[0]
        self.aa[0] = -1 * pD["xScale"] * ac # make attractive for U = -xScale * c       
         

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
        expression += '''+ {bb} * (y - {YY})^4'''.format(**fmt)
        # xExpression / gradient 
        #expression += '''+ {aa} * (x - {XX})^4'''.format(**fmt)
        # my potential for the x direction
        # xPotential = aa*(xParticle-x0)      <--- if bba=1
        expression += '''+ {aa} * (x - {XX})^{bba}  '''.format(**fmt)
        expression +=" "
  
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


##############################################################################
# Global parameters
##############################################################################

class Params():
  def __init__(self):
    paramDict = dict()

    #paramDict["friction"] = ( 50 / picosecond ) # rescaling to match exptl data PKH  
    paramDict["friction"] = ( 50              ) # rescaling to match exptl data PKH  
    paramDict["timestep"] = 10.0 * femtosecond# 1e-11 s --> * 100 --> 1e-9 [ns] 
                                #    72 s --> &*100 --> 7200 [s] (2hrs)     
    paramDict["nUpdates"] = 1000  # number of cycldes 
    paramDict["xPotential"] = False
    paramDict["yPotential"] = False
    paramDict["containmentPotential"] = False # 'cylindrical','square'
    paramDict["xScale"]   = 100.   # scale for chemoattractant gradient 
    paramDict["frameRate"]=   1.  # [1 min/update]

    paramDict["nCells"] = 10  
    paramDict["cellRad"] = 10.    # [um] cell radius  (seems too small) 
    paramDict["cellAttr"]=1.   # [] attraction between crowder and cell (vdw representation) 
    paramDict["cellChg"]= 0.   # [] 'electrostatic charge' (just for debugging)                        

    paramDict["nCrowders"] = 1  
    paramDict["crowderRad"]= 15. # [um]
    paramDict["crowderAttr"]=1.   # [] attraction between crowder and cell (vdw representation) 
    paramDict["crowderChg"]= 0.   # [] 'electrostatic charge' (just for debugging)                        
    paramDict["tag"]="run"
    paramDict["outName"]="test"

    # system params (can probably leave these alone in most cases
    paramDict["chemoAttr"]      = 1.     # conc. of chemoattractant 
    paramDict["domainDim"]    = 10  # FOR NOW, KEEP PARTICLES WITHIN 99 for PDB [nm/um] dimensions of domain  
    paramDict["crowderDim"]    = None   # [nm/um] dimensions of domain containing crowders (square)  
    paramDict["nInteg"] = 100  # integration step per cycle
    paramDict["mass"] = 1.0 * dalton
    paramDict["temperature"] = 298 * kelvin
    paramDict['dumpSize'] = 100 # make larger?  

    # store all values 
    self.paramDict = paramDict
     
  def update(self):
    if self.paramDict["crowderDim"] is None:
        self.paramDict["crowderDim"]=self.paramDict["domainDim"]

# allocate instance 
params = Params()

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


  crowderPos, cellPos = lattice.GenerateCrowdedLattice(
          nCrowders,nCells,
          paramDict['crowderRad'],paramDict['cellRad'],
          crowdedDim=paramDict["crowderDim"], # [um] dimensions of domain containing crowders (square)  
          outerDim=paramDict["domainDim"]
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
  #
  # START ITERATOR 
  #
  print("Running dynamics") 
  for i in range(nUpdates): 
      if (i % iters)==0:
        print("(%d/10)..."%(i/iters+1))

      # get positions at ith cycle 
      x = simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(nanometer)

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
  if diffMax > np.min([paramDict["crowderRad"],paramDict["cellRad"]]):
      print("Diff max %f"%diffMax)
      print("WARNING: step size is greater than particle sizes, which may give bad results;increase friction") 
  stop = timeit.default_timer()
  print('Time(s): ', int(stop - start))

  #if display:
  #    plt.show() 

  # package data 

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




