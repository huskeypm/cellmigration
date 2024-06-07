import numpy as np 
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
def GetContacts(dists,idxs,thresh=2):
    """
    Determines when two particles or more particles are colliding
    """
    s=dists[:,idxs]
    #print("sub",s)
    cellMin = np.min(s,axis=1)
    #idxClose = np.where(cellMin <= thresh) # will have one min distance for each cell
    #print('cellMin',cellMin)
    #print('iscloe',idxClose)
    #idxClose=1
    close = np.array(cellMin <= thresh,int)
    #print(close)
    return close

def CalcRDFs(traj,solvAtomName,soluteAtomName):
    """
    Iterates over all solute atoms to compute rdf 
    """
    solutes = []
    for res in traj.top.residues:
      if res.name == soluteAtomName:
          solutes.append( res.index + 1 )  # zero indexed

    if len(solutes)<1:
        raise RuntimeError(soluteAtomName + " not found ") 

    # just taker first for now
    #print(solutes) 
    mask1='@%s'%solvAtomName
    mask2=':%d@%s'%(solutes[0],soluteAtomName)

    CalcRDF(traj,mask1=mask1,mask2=mask2)

    return 

def CalcRDF(traj,
            mask1, # solvent
            mask2, # solute
            space=0.1,
            bins=20,
            display=False):
    import pytraj as pt
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
    bins,rdf = pt.rdf(traj, solvent_mask=mask1, solute_mask=mask2, 
            bin_spacing=1,maximum=500 # space, maximum=bins
            ) 

    #print(np.shape(radial))
    #print(radial)
    maxBin = bins[np.argmax(rdf)]

   
    if display:
      s = np.sum(rdf)
      plt.plot(bins*AA_to_NM,rdf/s)
      plt.xlabel("r [nm]") 
      plt.ylabel("P") 
      plt.title("RDF " + mask1 + " " + mask2)
      plt.gcf().savefig("rdf.png",dpi=300) 

    return maxBin 

def CalcProbDist(
  traj, mask='@RC',display=False,caseName=None,
    bins=None,   # binds in pixels (usually from px/um conversion) 
    tMax=-1):
  """
  Computes probability distribution and potential of mean force 
  (via boltzmann inversion) 
  for each particle, get dx in all directions, provide dt as input
  select particles in some neighborhood of y=0?
  """
  import pytraj as pt
  numFrames = (np.shape(traj.xyz))[0]

  ## get cells
  indices = pt.select_atoms(traj.top, mask)

  xs = traj.xyz[0:tMax,indices,0]
  xs = np.ndarray.flatten(xs)  # n particles x m timesteps 
  ys = traj.xyz[0:tMax,indices,1]
  ys = np.ndarray.flatten(ys)  # n particles x m timesteps 

  if bins is not None:
    p,x,y= np.histogram2d(xs,ys,bins=bins,density=True)
  else:
    p,x,y= np.histogram2d(xs,ys,density=True)
  np.savetxt(caseName+"prob.csv",p) 

  dx=x[1]-x[0]
  dy=y[1]-y[0]
  X, Y = np.meshgrid(x, y)

  # get PMF
  p[p<thresh]=thresh
  pmf = -np.log(p) * kT 
  pmf-=np.min(pmf) 


  if caseName is None:
    caseName=""
  else: 
    caseName+="_"

  #display=True
  if display:
    plt.figure()
    plt.axis('equal')
    plt.pcolormesh(X, Y, p.T) # probably T is appropriate here 
    plt.colorbar()
    plt.gcf().savefig(caseName+"prob2d.png",dpi=300)

    plt.figure()
    plt.axis('equal')
    plt.pcolormesh(X, Y, pmf.T) # probably T is appropriate here 
    plt.colorbar()
    plt.gcf().savefig(caseName+"pmf.png",dpi=300)

  return p,X,Y,dx,dy

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
def LoadTraj(caseName,warningOnly=False):
    import pytraj as pt 
    # load
    try:
      dcd=caseName+".dcd"; pdb=caseName+".pdb"
      traj = pt.iterload(dcd,pdb)
    except:
      if warningOnly:
        print("You're likely missing a file like %s"%dcd)
        return None 
      else: 
        raise RuntimeError("You're likely missing a file like %s"%dcd)
    print("Loaded %s"%dcd)

    return traj

def CalcVolFrac(auxParams): 
  """ 
  assumes vol frac is defined by domainYDim (domainYDim = crowderDim + crowderRad)
  """

  d = auxParams['domainYDim']
  n = auxParams['nCrowders']

  r = auxParams['crowderRad']
  if('effectiveRad' in auxParams.keys() and auxParams['effectiveRad'] is not None):
      print("using eff. radius for vol frac") 
      r = auxParams['effectiveRad']

  volFrac = CalcVolFraction(n,d,r)

  return volFrac 

def CalcVolFraction(
    nCrowder,
    dimRegion=200.,
    crowderRad=10.):
    """
    Calculates the volume fraction occupied by crowders
    Works better once an effective radius is computed via simulatons/probability distribution
    """
    
    areaCrowder= nCrowder *np.pi * crowderRad **2
    areaRegion= dimRegion**2
    volFrac= (areaRegion-areaCrowder)/areaRegion
    return volFrac

def CalcD(traj,mask='@RC',csvName=None, display=False):
  """
  Calculates diffusion coefficient via mean squared displacements
  """
  import pytraj as pt 
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
    plt.plot(ts_MIN[equilFrame:],msd_NMNM[equilFrame:])
    plt.gcf().savefig(csvName+".png")

  # x**2 = 4*D*t
  #      = slope * t ==> D = slope/4.
  Di = slope/4.
  return Di 

# assumes cylindrical inclusions 

# J = D * grad(c)
# derivative in x direction 
#diff = np.abs( np.diff(data2,axis=0) ) # to check that zeroing out works 
from scipy.ndimage import gaussian_filter
from scipy.stats import linregress
def CalcAverageFlux(
    dataSet,
    D=1,
    barrierLims = [50,60],
    reservoirLims = [50,60],
    dx = 0.125,
    dy=  0.125,
    display=False):# need to relate px to nm
    """
    Calculates flux over region in order to get average values. Considers both crowder location and left reservoir
    dataSet -  # concentration/probability
    D - diffusion coefficient
    barrierLims - range in x over which to evaluate D in crowded region
    reservoirLims - range in x over which to evaluate D in reservoir region
    dx,dy - resolution in nm per px (currently 1400 px wide for 350 nm domain)
    """
    #data2 = np.outer(-np.arange(100),np.ones(100))
    #plt.pcolormesh(data2.T)

    # make mask
    mask = np.ones_like(dataSet)
    avgDataSet = np.mean(dataSet)
    dataSetCp = np.copy(dataSet)
    #print('data set avg', avgDataSet)
    cutoff = 0.3*avgDataSet
    mask[dataSetCp < cutoff ]= 0   # 
    insideMean = np.mean(dataSetCp[dataSetCp < cutoff ]) # mean inside inclusion 
    outsideMean = np.mean(dataSetCp[dataSetCp >= cutoff ]) # mean outside inclusion
    
    cutoff = 0.1*avgDataSet
    dataSetCp[dataSetCp < cutoff]= 0 # outsideMean   # we apply a more stringent criterion here so that we don't 'zero' out as much of the crowder-->smoother gradients
    
    mask=mask[1:,:] #trim first row to match with diff later
    
    #plt.pcolormesh(mask.T) # may be zero if no occlusions 
    if display:
      plt.figure()
      fig, axs = plt.subplots(1)
      axs.pcolormesh(mask.T) #. mask[xlims[0]:xlims[1]].T)
      axs.set_aspect('equal')       
      axs.set_title("mask")
    
    # smoothing 
    
    #dataSetCp = gaussian_filter(dataSetCp, sigma=3)
    Dgradc = D*np.diff(dataSetCp,axis=0)  # do verify yhis is the right direction 
 
    # make out low sampled area (occlusions)
    #print(np.mean(diff))
    J=Dgradc*mask
    #plt.hist(J)

    #display=True
    xlims = barrierLims 
    if display:
      plt.figure()
      fig, axs = plt.subplots(1)
      axs.pcolormesh(dataSetCp[xlims[0]:xlims[1]].T)
      axs.set_aspect('equal')       
      axs.set_title("modified c")

      plt.figure()
      fig, axs = plt.subplots(1)
      axs.pcolormesh(J[xlims[0]:xlims[1]].T)
      axs.set_aspect('equal')       
      axs.set_title("masked J")
      plt.gcf().savefig("avgflux.png",dpi=1200)
    

    # get subregion containing the occlusions 
    subJ =J[xlims[0]:xlims[1],:]
    submask =mask[xlims[0]:xlims[1],:]

    nx,ny = np.shape(submask)
    #plt.pcolormesh((submask).T)
    #plt.pcolormesh((subJ).T)
    JSum = np.sum(subJ)*dx*dy  # J = D * grad(c)
    areaSum = np.sum(submask)*dx*dy
    areaTot = (dx*nx)*(dy*ny)
    #print(nx*ny)
    JavgCrowded = JSum/areaTot
    areaFrac = areaSum/areaTot
    #print(areaFrac,Javg)
    
    # try part before crowders too
    subset=dataSetCp #[:,50:350]
    dasum = np.sum(subset,axis=1)
    dasum = gaussian_filter(dasum, sigma=3)    
    if display:
        plt.figure()
        fig, axs = plt.subplots(1)
        #subset=Dgradc   # [100:500,50:350]
        axs.pcolormesh(subset.T,cmap='gray')
    
    # look in first 1/4 of plot to find max     
    #upper = int( (xlims[0])/3.) # if slope is -increasing- [happens w ATP], we can get the false limit
    #daMax = np.argmax(dasum[0:upper])
    #lims=[daMax,xlims[0]]
    lims = reservoirLims
    x = np.linspace(lims[0],lims[1],lims[1]-lims[0])
    y = dasum[lims[0]:lims[1]]
    result = linregress(x, y)
    JavgReservoir=result.slope
    
    #print("JavgR ",JavgReservoir," JavgCrowd ", JavgCrowded, " areafrac ", areaFrac)
    
    #display=True
    if display:
        plt.figure()
        plt.plot(dasum)
        plt.plot(x, x*result.slope + result.intercept,label=result.slope)
        plt.legend(loc=0)
        plt.gcf().savefig("testsum.png",dpi=300)
    
    
    return areaFrac,JavgCrowded, JavgReservoir




## get flux
def CalcFluxLine(traj, mask='@RC',display=False,xThresh=0,margin=None,caseName=None): 
  """
  Gets particle flux across xThresh (0 for now) 
  for each particle, get dx in all directions, provide dt as input
  select particles in some marging of xThresh   
  """

  ## get cells 
  import pytraj as pt
  indices = pt.select_atoms(traj.top, mask) 
  
  # xyz: frames, natoms, coords
  #print(traj.xyz[2,:,0])

  # i can spread this out over an interval
  if display: 
    xdiffs = traj.xyz[-1,indices,0] 
    xdiffs -= traj.xyz[0,indices,0] 
    plt.figure()
    #print(traj.xyz[0,indices,0])
    #print(traj.xyz[-1,indices,0])
    plt.xlabel('P(xdisplacements)')
    plt.title("XDisplacements(tf-t0)") 
    plt.hist(xdiffs)
    plt.gcf().savefig("diffs.png") 
    print("Mean displacement",np.mean(xdiffs))
  


  
  ## below Thresh (only need to track one 'compartment' for now since we are dividing the space into two sections
  ## along x direction 
  ##l =np.array(traj.xyz[:,0,0] < xThresh, dtype=int) 
  ##r =np.array(traj.xyz[:,0,0] >=xThresh, dtype=int) 
  ##print(np.sum(l),np.sum(r))
  ## indicates which of the particles are left of xThresh
  #l =np.array(traj.xyz[:,indices,0] < xThresh, dtype=int) 
  ## counts number of particles left of xThresh for each time step 
  #l =np.sum(l,axis=1)
  ## if the population in the compartment changes between two times, then a particle has entered/left
  #diff = np.diff(l)  
  #fluxArea = diff/dt # not normalizing by area
  #JA = np.average(fluxArea)

  
  # for testing 
  if False: 
    n = 5
    z = np.outer(np.arange(n) + 0.01*np.random.randn(n),[1,1]) 
    z[:,1] = n - np.arange(n) + 0.01*np.random.randn(n)
    lz=np.array(
    [[1,5,0,3,5],    
     [2,4,2,2,5], 
     [3,3,4,3,2],
     [4,2,2,2,2], 
     [5,1,0,3,2]]
    )
    #print(np.shape(z))
  xs = traj.xyz[:,indices,0]      
  zshifted = xs - xThresh  # values to the right of xThresh are positive; otherwise negative
  print("Flux using xThresh ",xThresh)
  
  # exclude those outside of margin 
  if margin is not None: 
    zshifted[ np.abs(zshifted)> margin] = 0
    #print(zshifted)
  
  
  p = zshifted[1:,] * zshifted[0:-1,] # if the product is negative, the particle has crossed the threshold (l or r)
  #print(p)
  #negIdx = np.argwhere(p<0)
  #print(negIdx) # just for easy visualization 
  negIdx = np.where(p<0)
  changed = np.zeros_like(p)
  changed[negIdx]=1
  #print(changed)
  
  zshifted = zshifted[1:,] 
  logic = zshifted * changed # if negative (e.g. zshift<0, p=1),particle moved left
                                  # else particle moved right
  changed[ logic < 0 ] *= -1
  
  
  #print(changed)
  netCrossed=np.sum(changed,axis=1)
  #print(netCrossed)
  
  


  # get cumulative mean 
  timeSoFar = np.arange( np.shape(netCrossed)[0] )+1
  cummean = np.cumsum(netCrossed)/timeSoFar
  JA = cummean[-1]
  
  l =np.array(traj.xyz[:,indices,0] < xThresh, dtype=int) 
  l =np.sum(l,axis=1)
  if (l[-1] < 0.1*l[0]):
      print("WARNING: compartment is nearly depleted/flux estimates may be unreliable") 
  #print("J*A = %f "%JA)
  
  #x = traj.xyz[:,indices,0]
  #y = traj.xyz[:,indices,1]
  #plt.plot(x,y)
  #plt.gcf().savefig("testxy.png") 
  if caseName is None:
    caseName=""
  if display:
    plt.figure()
    print("Flux %8f"%cummean[-1] )
    axl = plt.subplot(111)
    #axl.plot(fluxArea,label="flux*area")
    axl.plot(timeSoFar[100:],cummean[100:],label="flux*area")
    #axl.set_xlim(0,2500)
    axr = axl.twinx()
    axr.plot(l,'r',label="#particles in x<thresh")     
    axr.set_ylim(0,np.max(l)+1)
    axl.legend(loc=1)
    axr.legend(loc=2)
    plt.gcf().savefig(caseName+"flux.png",dpi=600) 

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
        

