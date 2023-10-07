import matplotlib.pylab as plt
import numpy as np 
from sklearn.linear_model import LinearRegression

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

def CalcVolFrac(
    nCrowder,
    dimRegion=200.,
    crowderRad=10.):
    
    areaCrowder= nCrowder *np.pi * crowderRad **2
    areaRegion= dimRegion**2
    volFrac= (areaRegion-areaCrowder)/areaRegion
    return volFrac

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
