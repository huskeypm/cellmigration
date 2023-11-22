import numpy as np 
import pytraj as pt
import matplotlib.pylab as plt

equilFrame = 400
equilFrame = 0
FRAME_CONV = 0.1 # min/fr 
dt = FRAME_CONV  # [min] 
kT = 1  # [kT] 
thresh=1e-9
AA_to_NM = 0.1

## 
## FUNC
##
def CalcRDF(traj,
            mask1, # solvent
            mask2, # solute
            space=0.1,
            bins=20):
    #radial = pt.rdf(traj, solvent_mask=':WAT@O', solute_mask=':WAT@O', bin_spacing=0.2, maximum=12.)
    """ 
    Calculates radial distribution function
    traj
    mask1 - solvent
    mask2 - solute 
    """
    # works 
    #mask1='@1-20'    
    #mask2='@28'
    # works 
    #mask1='@RC'    
    #mask2=':28@AC'

    # rdf 
    radial = pt.rdf(traj, solvent_mask=mask1, solute_mask=mask2, 
            bin_spacing=1,maximum=500 # space, maximum=bins
            ) 

    #print(np.shape(radial))
    #print(radial)

   
    s = np.sum(radial[1])
    plt.plot(radial[0]*AA_to_NM,radial[1]/s)
    plt.xlabel("r [nm]") 
    plt.ylabel("P") 
    plt.title("RDF " + mask1 + " " + mask2)
    plt.gcf().savefig("rdf.png",dpi=300) 
    return radial 

def CalcProbDist(traj, mask='@RC',display=False):
# for each particle, get dx in all directions, provide dt as input
# select particles in some neighborhood of y=0?
# grab dx along flux direction
# sum(dx) / delta y
  numFrames = (np.shape(traj.xyz))[0]

  ## get cells
  indices = pt.select_atoms(traj.top, mask)

  xs = traj.xyz[0:,indices,0]
  xs = np.ndarray.flatten(xs)  # n particles x m timesteps 
  ys = traj.xyz[0:,indices,1]
  ys = np.ndarray.flatten(ys)  # n particles x m timesteps 

  bins = 100
  p,x,y= np.histogram2d(xs,ys,bins=bins,density=True)
  X, Y = np.meshgrid(x, y)

  # get PMF
  p[p<thresh]=thresh
  pmf = -np.log(p) * kT 
  pmf-=np.min(pmf) 


  #display=True
  if display:
    plt.figure()
    plt.axis('equal')
    plt.pcolormesh(X, Y, p.T) # probably T is appropriate here 
    plt.colorbar()
    plt.gcf().savefig("prob.png",dpi=300)

    plt.figure()
    plt.axis('equal')
    plt.pcolormesh(X, Y, pmf.T) # probably T is appropriate here 
    plt.colorbar()
    plt.gcf().savefig("pmf.png",dpi=300)

  return p,X,Y

def PlotStuff(
  msds,
  ts,
  outName=None
  ):     
  raise RuntimeError("PLOTSTUFF has been renamed; use CalcMSD instead\n") 
  CalcMSD(ts,msds,outName=outName)

def meanSquareDisplacements(xs, ys, nUpdates):
    x0s = xs[:,0]
    y0s = ys[:,1]
    msds = np.zeros(nUpdates)
    for i in range(nUpdates):
      xmsd = xs[:,i] - x0s
      ymsd = ys[:,i] - y0s

      sd = xmsd*xmsd + ymsd*ymsd
      msd = np.mean(sd)
      msds[i] = msd

    return msds

# casename should include full path 
def LoadTraj(caseName):
    # load
    try:
      dcd=caseName+".dcd"; pdb=caseName+".pdb"
      traj = pt.iterload(dcd,pdb)
    except:
      raise RuntimeError("You're likely missing a file like %s"%dcd)
    print("Loaded %s"%dcd)
    return traj

def CalcVolFrac(auxParams): 
  d = auxParams['domainDim']
  n = auxParams['nCrowders']

  if auxParams['effectiveRad'] is None:
    r = auxParams['crowderRad']
  else:
    r = auxParams['effectiveRad']
    print("USING EFFECTIVE RADIUS") 

  volFrac = CalcVolFraction(n,d,r)

  return volFrac 

def CalcVolFraction(
    nCrowder,
    dimRegion=200.,
    crowderRad=10.):
    
    areaCrowder= nCrowder *np.pi * crowderRad **2
    areaRegion= dimRegion**2
    volFrac= (areaRegion-areaCrowder)/areaRegion
    return volFrac

def CalcD(traj,mask='@RC',csvName=None, display=False):
  # in A^2 
  rmsd = pt.rmsd(traj, mask='@RC', ref=0)
  #print(rmsd[0:10])
  msd = rmsd**2
  # in NM^2
  AA_to_NMNM=1e-2
  msd_NMNM = msd*AA_to_NMNM

  # in XXX min/fr
  tEnd = np.shape(rmsd)[0]
  ts_MIN = np.arange(tEnd) * dt

  # fit for D [nm^2/min]
  slope,intercept= np.polyfit(ts_MIN[equilFrame:],msd_NMNM[equilFrame:],1)
  #print(slope)
  #plt.plot(rmsd)
  #plt.plot(ts,ts*slope+intercept)
  #

  if csvName is not None:
    dim = np.shape(rmsd)[0]
    csv = np.reshape(np.zeros(dim*2),[dim,2])
    csv[:,0] = ts_MIN 
    csv[:,1] = msd_NMNM
    np.savetxt(csvName+".csv",csv,header='time[min],msd[nmnm]')   

  if display: 
    plt.figure()
    plt.plot(ts,rmsd)
    plt.gcf().savefig(csvName+".png")

  # x**2 = 4*D*t
  #      = slope * t ==> D = slope/4.
  Di = slope/4.
  return Di 

# assumes cylindrical inclusions 


## get flux
def CalcFlux(traj, mask='@RC',display=False):
# for each particle, get dx in all directions, provide dt as input
# select particles in some neighborhood of y=0?
# grab dx along flux direction
# sum(dx) / delta y
  numFrames = (np.shape(traj.xyz))[0]

  ## get cells 
  indices = pt.select_atoms(traj.top, mask) 
  
  # xyz: frames, natoms, coords
  #print(traj.xyz[2,:,0])
  xThresh = -0.

  if display: 
    plt.figure()
    diffs = traj.xyz[-1,indices,0] 
    diffs -= traj.xyz[0,indices,0] 
    #print(traj.xyz[0,indices,0])
    #print(traj.xyz[-1,indices,0])
    plt.hist(diffs)
    plt.gcf().savefig("diffs.png") 
  


  
  # below Thresh (only need to track one 'compartment' for now since we are dividing the space into two sections
  # along x direction 
  #l =np.array(traj.xyz[:,0,0] < xThresh, dtype=int) 
  #r =np.array(traj.xyz[:,0,0] >=xThresh, dtype=int) 
  #print(np.sum(l),np.sum(r))
  l =np.array(traj.xyz[:,indices,0] < xThresh, dtype=int) 
  l =np.sum(l,axis=1)
  # if the population in the compartment changes between two times, then a particle has entered/left
  diff = np.diff(l)  
  fluxArea = diff/dt # not normalizing by area
  JA = np.average(fluxArea)
  if (l[-1] < 0.1*l[0]):
      print("WARNING: compartment is nearly depleted/flux estimates may be unreliable") 
  #print("J*A = %f "%JA)
  
  #x = traj.xyz[:,indices,0]
  #y = traj.xyz[:,indices,1]
  #plt.plot(x,y)
  #plt.gcf().savefig("testxy.png") 
  if display:
    print(l[0],l[-1],np.sum(fluxArea),np.average(fluxArea)*numFrames) # 320 frames 
    plt.plot(l,label="#particles in x<thresh")     
    plt.plot(fluxArea,label="flux*area")
    plt.legend(loc=0)
    plt.gcf().savefig("test.png") 

  return JA         





def CalcMSD(     
  ts,
  msds,
  outName=None,
  display=True
  ):     
  """
  ts - len(ts) = num frames from microscope (assuming 1/min) 
  msds - [um]  
  """
  raise RuntimeError("phase out") 

  from sklearn.linear_model import LinearRegression
  #print(np.shape(ts))
  #print(ts[-1])

  # display rmsd
  # exp data later 
  adjTime = 1.
  texp = ts*adjTime
  texp_=texp.reshape((-1, 1))

  model = LinearRegression().fit(texp_,msds)
  #model = LinearRegression().fit(texp,msds)
  msdfit= model.intercept_ + texp*model.coef_
  D = model.coef_ / 4.   # since msd = 4Dt for 2D diffusion 

  import matplotlib.pylab as plt
  if display:
    plt.plot(texp,msds  , label="exp")
    plt.plot(texp,msdfit, label="fit")
    plt.ylabel("MSD [um^2]")
    plt.xlabel("time [min]")
    plt.title("D= %f [um^2/min]"%D) 
    plt.legend(loc=0)
 
    # 1K ts --> 2hrs/7200s?  
    # 200 [um**2 ] 
    #plt.show()
    #plt.gca().savefig("compare.png")
    if outName is not None:
      plt.savefig(outName)       

  return texp,msdfit,D 


def PlotFinalPosition(xs,ys,outName=None):
    sx=xs[:,0] 
    sy=ys[:,0]
    plt.scatter(sx,sy,label='start')
    fx=xs[:,-1]
    fy=ys[:,-1]
    plt.scatter(fx,fy,label='final')

    coms=np.array([np.mean(sx),np.mean(sy)])
    comf=np.array([np.mean(fx),np.mean(fy)])
    #print(coms, comf)
    dist = (comf - coms)**2
    dist = np.sqrt( np.sum(dist) )
    plt.title("COM displacement %f"%dist)


    if outName is None:
        outName = "out.png"
    plt.gcf().savefig(outName)

    


import pickle as pkl

def LoadPKLData(trajOutName="test.pkl"):
    file = open(trajOutName,"rb")
    data = pkl.load(file)
    file.close()

    [ts,xs,ys] = data
    nUpdates = np.shape(ts)[0]
    nParticles = np.shape(xs)[0]

    return ts,xs,ys,nUpdates,nParticles
###                           PDB Generator                             ###
###########################################################################

def genPDBWrapper(
        pdbFileName,
        nCellsType1,   # generally the diffusing cells 
        nCellsType2=0, # generally the crowders 
        startingPositions=None):
    """
    This is a wrapper for Ben's PDB writer, siunce I needed something simple.
    Later I can revise his code to accept coordinates; currently they're randomized 
    """
    #print(startingPositions) 
    PDBgenNoPBC(pdbFileName,nCellsType1,nCellsType2,0,0,"None",startingPositions=startingPositions)
    

def PDBgenNoPBC(PDBfileName,
                NoCell1, 
                NoCell2,
                numOfDeadCell,
                num_BP,
                DiffState,  
                startingPositions=None):
    '''
    
    [Note]: 11/06/2020 by Ben 
    In this function, having all resting cells or activated cells will not cause any error.
    
    [Parameter Description]
    PDBfileName:    pdb file name
    NoCell1:        a total number of cells in the 1st half of box: should be resting cells
    NoCell2:        a total number of cells in the 2nd half of box: should be activated cells 

    '''    
    # --- Writing pdb file initiated ---
    daFile = PDBfileName
    if ".pdb" not in daFile:
        daFile = daFile+".pdb"
    structure = open(daFile,"w") 
    structure.write("MODEL     1 \n")

    '''
    [------First half of box-------]
    '''
    # determining a total number of resting cells from the overall population
    ## the rest of cells are activated cells
    
    TotalNoCell = NoCell1 + NoCell2 + numOfDeadCell + num_BP
    refnum1 = NoCell1 + NoCell2
    refnum2 = refnum1 + num_BP
    refnum3 = numOfDeadCell + NoCell1
    refnum4 = refnum3 + NoCell2
    
    ## Dead Cell first: There is no distinction between 1st and 2nd compartment for the dead cells 
    for i in np.arange(TotalNoCell):
        # add random positions if starting positions are not defined 
        if startingPositions is None:
          x = randint(0,9)
          y = randint(0,9)
          z = randint(0,9)
        else:
          #print(np.shape(startingPositions))
          x,y,z = startingPositions[i,:]

        # make sure num sig. figs is correct for pdb
        x = format(x,'8.3f') # any changest to this need to be reflected in spacing below 
        y = format(y,'8.3f')
        z = format(z,'8.3f')
        
        if numOfDeadCell == 0:
            if i < NoCell1:
                name = 'RC'
            elif NoCell1 <= i < refnum1 :
                if DiffState =='steady':
                    name = 'RC'
                else:
                    name = 'AC'
            else:
                name = 'BC'
        else:
            if i < numOfDeadCell:
                name = 'DC'
            elif numOfDeadCell <= i < refnum3: 
                name = 'RC'
            elif refnum3 <= i < refnum4:
                if DiffState == 'steady':
                    name = 'RC'
                else:
                    name = 'AC'
            else:
                name = 'BC'

        if i < 9:
            structure.write("ATOM      "+str(int(i+1))+"  "+name+"   "+name+"     "    +str(int(i+1))+"    "+str(x)+""+str(y)+""+str(z)+"  1.00  0.00 \n")
        elif i >= 9 and i < 99:
            structure.write("ATOM     "+ str(int(i+1))+"  "+name+"   "+name+  "    "   +str(int(i+1))+"    "+str(x)+""+str(y)+""+str(z)+"  1.00  0.00 \n")
        elif i >= 99 and i < 999:
            structure.write("ATOM    "+  str(int(i+1))+"  "+name+"   "+name+    "   "  +str(int(i+1))+"    "+str(x)+""+str(y)+""+str(z)+"  1.00  0.00 \n")
        elif i >= 999 and i < 9999:
            structure.write("ATOM   "+   str(int(i+1))+"  "+name+"   "+name+      "  " +str(int(i+1))+"    "+str(x)+""+str(y)+""+str(z)+"  1.00  0.00 \n")
        elif i >= 9999:
            structure.write("ATOM  "+    str(int(i+1))+"  "+name+"   "+name+        " "+str(int(i+1))+"    "+str(x)+""+str(y)+""+str(z)+"  1.00  0.00 \n")
            
    structure.write("ENDMDL")
    structure.close
                   
    return
        

