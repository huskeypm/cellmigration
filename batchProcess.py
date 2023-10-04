#!/usr/bin/env python
import matplotlib.pylab as plt
import numpy as np
import brown_util as bu

# pip3 install pytraj
import pytraj as pt


# In[5]:


# insert and load the names of the pkl files that you made using the yamlFile to change the parameters
# calculate the mean square displacements for each pkl file
# plot the simulated trajectories of the particles


cases=dict()
trajNames = ["10","20"]
#use dcd
equilFrame = 400
iters=1

slopes = []
for trajName in trajNames:
  for iter in range(iters):

    caseName = trajName + "_%d"%(iter)
    print(caseName)
    next
  
    # load
    traj = pt.iterload(caseName+".dcd", caseName+".pdb")

    # get rmsd 
    mask='@RC'
    rmsdAll = pt.rmsd(traj, mask='@RC', ref=0)
    rmsd = rmsdAll[equilFrame:]
    tEnd = np.shape(rmsd)[0] 
    dt = 1.
    ts = np.arange(tEnd) * dt
    
    # fit for D 
    #slope,intercept= np.polyfit(ts, ts*iter, 1)
    slope,intercept= np.polyfit(ts, rmsd, 1)
    slopes.append(slope)
    #print(slope)
    #plt.plot(rmsd)
    #plt.plot(ts,ts*slope+intercept)
    # 

  # x**2 = 4*D*t
  #      = slope * t ==> D = slope/4.
  Ds =  np.asarray( slopes )/4.

  case = empty()
  case.Dmean = np.mean( Ds )
  case.Dstd = np.std( Ds)

  cases[pathRatio] = case
quit()


# In[51]:



Ds = []
DsErr = []
inds = range( len( pathRatios ) )
for pathRatio in pathRatios:
  case = cases[pathRatio]
  Ds.append(case.Dmean)
  DsErr.append(case.Dstd)
  
# normalize by fastest Ds (widest)
nDs = np.array( Ds ) / Ds[-1]
nPRs = np.array( pathRatios ) / pathRatios[-1]

print(Ds)
#plt.bar(inds,Ds) 
#print(nDs)
plt.plot(nPRs,nDs,'k.',label="Sims")
plt.plot(nPRs,nDs,'k')
plt.plot(nPRs,nPRs,'r.',label="Analytic")


# In[ ]:





quit()
path="tests/"
for trajectory in listofTrajectories:
    ts,xs,ys, nUpdates, nParticles = bu.LoadPKLData(path+trajectory)
    msds= bu.meanSquareDisplacements(xs, ys, nUpdates)
    #print(msds[-1])
    #plotTrajectories(xs, ys)
    #plt.gcf().savefig("testing_NumberofCrowders_particles.png") #how to save each of these individually?


# In[7]:


allData=[]
listNames= ["crowder_EK_16_crowders", "crowder_EK_25_crowders"]
listnCrowders=[16,25]
nRuns= 10
path="/home/ekrueger2/source/cellmigration/tests/"

for i, name in enumerate(listNames):
    ds=np.zeros(nRuns)
    for j in range(nRuns):
        namei= name + "_%d"%(j+1) + ".pkl"
        print(namei)
        ts,xs,ys, nUpdates, nParticles = bu.LoadPKLData(path+namei)
        msds= bu.meanSquareDisplacements(xs, ys, nUpdates)
        #bu.PlotStuff(msds,ts)
        texp, msdexp, Dexp = bu.CalcMSD(ts,msds)
        print(Dexp)
        ds[j]=Dexp
    datai={'Dexp_mean':np.mean(ds), 'Dexp_sd':np.std(ds)}
    allData.append(datai)
    plt.gcf().savefig("16_and_25_crowders_10runs.png")


# In[9]:


ds=np.zeros(len(datai))
sds=np.zeros(len(datai))

for i, datai in enumerate(allData):
    ds[i]=datai['Dexp_mean']
    sds[i]=datai['Dexp_sd']
#inds=np.arange(len(datai))
inds=bu.CalcVolFrac(np.asarray(listnCrowders))
import matplotlib.pylab as plt
plt.scatter(inds, ds)
plt.errorbar(inds, ds, yerr=sds)
plt.gcf().savefig("16_and_25_crowders_10runs_errorBars.png")


# In[ ]:




