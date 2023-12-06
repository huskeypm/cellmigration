##############################################################################
# Global parameters
##############################################################################
from simtk.unit import kelvin, picosecond, femtosecond, nanometer, dalton

class Params():
  """
  Class the contains paramDict object with program parameters"
  """
  def __init__(self):
    paramDict = dict()
    #paramDict["friction"] = ( 50 / picosecond ) # rescaling to match exptl data PKH  
    paramDict["friction"] = ( 50              ) # rescaling to match exptl data PKH  
    paramDict["timestep"] = 10.0 * femtosecond# 1e-11 s --> * 100 --> 1e-9 [ns] 
                                #    72 s --> &*100 --> 7200 [s] (2hrs)     
    paramDict["nUpdates"] = 1000  # number of cycldes 
    paramDict["xPotential"] = False # enable potential in x direction 
    paramDict["yPotential"] = False
    paramDict["containmentPotential"] = False # 'cylindrical','square'
    paramDict["chemoAttr"]      = 1.     # conc. of chemoattractant 
    paramDict["xScale"]   = 10.   # scale for chemoattractant gradient along X 
    paramDict["yScale"]   = 10.   # scale for chemoattractant gradient along Y 
    paramDict["frameRate"]=   1.  # [1 min/update]

    paramDict["nCells"] = 10  # number of cells  
    paramDict["cellRad"] = 10.    # [um] cell radius  (seems too small) 
    paramDict["cellAttr"]=1.   # [] attraction between crowder and cell (vdw representation) 
    paramDict["cellChg"]= 0.   # [] 'electrostatic charge' (just for debugging)                        

    paramDict["nCrowders"] = 1     # number of crowders 
    paramDict["crowderRad"]= 15. # radius of crowders [um]
    paramDict["crowderAttr"]=1.   # [] attraction between crowder and cell (vdw representation) 
    paramDict["crowderChg"]= 0.   # [] 'electrostatic charge' (just for debugging)                        
    paramDict["contactDistCrowderA"] = 1000.  # distance between cell and crowder type A to be considered a 'contact'
    paramDict["contactDistCrowderB"] = 1000. # distance between cell and crowder type A to be considered a 'contact'
    paramDict["effectiveRad"] = None    # [nm] this is used when the attraction between crowder/cell yields effective radii that are smaller than expected.  
    paramDict["tag"]="run" # optional tag name for job 
    paramDict["outName"]="test" # name for ouput files 

    # for states
    paramDict["states"] = False  # turn on states calculations(see states.py for more info)   
    paramDict["K01"] = 1000.  # rate for 0->1 transition in states.py (see that file for more info)   
    paramDict["K12"] = 1000.  # rate for 0->1 transition in states.py (see that file for more info)   
    paramDict["K20"] = 1000.  # rate for 0->1 transition in states.py (see that file for more info)   

    # system params (can probably leave these alone in most cases
    paramDict["domainXDim"]    = 50  # FOR NOW, KEEP PARTICLES WITHIN 99 for PDB [nm/um] dimensions of domain  
    paramDict["domainYDim"]    = 50  # FOR NOW, KEEP PARTICLES WITHIN 99 for PDB [nm/um] dimensions of domain  
    paramDict["xThresh"]       = None  # 0 keep only the particles on the left side 
    paramDict["absorbingBoundary"]     = False  # if True, when a particle is within absrbingMargin of the right boundary, it will be moved to the left boundary 
    paramDict["absorbingMargin"]       = 2
    paramDict["crowderXDim"]    = None   # [nm/um] dimensions of domain containing crowders (square)  
    paramDict["crowderYDim"]    = None   # [nm/um] dimensions of domain containing crowders (square)  
    paramDict["nInteg"] = 100  # integration step per cycle
    paramDict["mass"] = 1.0 * dalton  # particle mass 
    paramDict["temperature"] = 298 * kelvin # system temperature [K] 
    paramDict['dumpSize'] = 100 # make larger?  

    # store all values 
    self.paramDict = paramDict
     
  def update(self):
    if self.paramDict["crowderYDim"] is None:
        # set to Y-dim, which should usually be smallest
        self.paramDict["crowderYDim"]=self.paramDict["domainYDim"]
        self.paramDict["crowderXDim"]=self.paramDict["domainYDim"]


