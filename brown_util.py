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
def GetContacts(dists,idxs,thresh=2):
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
    Iteratures over all solute atoms to compute rdf 
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
  """
  Gets particle flux across xThresh (0 for now) 
  for each particle, get dx in all directions, provide dt as input
  select particles in some neighborhood of y=0?
  grab dx along flux direction
  sum(dx) / delta y
  """
  numFrames = (np.shape(traj.xyz))[0]

  ## get cells 
  indices = pt.select_atoms(traj.top, mask) 
  
  # xyz: frames, natoms, coords
  #print(traj.xyz[2,:,0])
  xThresh = -0.

  # i can spread this out over an interval
  xdiffs = traj.xyz[-1,indices,0] 
  xdiffs -= traj.xyz[0,indices,0] 
  if display: 
    plt.figure()
    #print(traj.xyz[0,indices,0])
    #print(traj.xyz[-1,indices,0])
    plt.xlabel('P(xdisplacements)')
    plt.title("XDisplacements(tf-t0)") 
    plt.hist(xdiffs)
    plt.gcf().savefig("diffs.png") 
  print("Mean displacement",np.mean(xdiffs))
  


  
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
    plt.figure()
    print(l[0],l[-1],np.sum(fluxArea),np.average(fluxArea)*numFrames) # 320 frames 
    axl = plt.subplot(111)
    axl.plot(fluxArea,label="flux*area")
    axr = axl.twinx()
    axr.plot(l,'r',label="#particles in x<thresh")     
    axr.set_ylim(0,np.max(l)+1)
    axr.legend(loc=0)
    plt.gcf().savefig("flux.png",dpi=600) 

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
