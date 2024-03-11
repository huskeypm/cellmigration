import matplotlib.pylab as plt
import numpy as np
import trackpy as tp 
from matplotlib import rcParams
import pandas as pd


def TrialParams(frames,
              refFrame=0, # reference frame to find paticles
              diameter=41, # diameter of the particle in pixel, must be an odd number
              minmass=1e2, # total brightness of a particle: a threshold to discard floating particles
             ):

    """
    Get tracking params from a given reference frame
    """
    rcParams['figure.figsize']=3,3
    # first pass
    f = tp.locate(frames[refFrame], diameter)
    fig, ax = plt.subplots()
    ax.set_title('first trial (diam only)')
    tp.annotate(f, frames[refFrame])

    # plot mass distribution
    fig, ax = plt.subplots()
    ax.set_title('mass histogram')
    ax.hist(f['mass'],bins=20)

    # filter using mass
    f = tp.locate(frames[refFrame], diameter, minmass=minmass)
    fig, ax = plt.subplots()
    ax.set_title('filtered by mass')
    tp.annotate(f, frames[refFrame])

    return f

def DoTracking(frames,diameter=41,minmass=1e2):
    """
    track particles using params determined from getParams()
    """
    fb = tp.batch(frames, diameter, minmass=minmass)
    return fb

def DoMSD(
    fb,
    maxDist=10, # maximum displacement between frames in pixel
    maxMissFrame=200, # allowed number of frames a particle can disappear
    minFrame=50, # minimum number of frames a trajectory needs to last
    name=None, # specify file name if to save trajectory figure
    pixelSize=1.63, # image pixel size in micron/pixel, specific for each scope and lens (Cytiva: 0.65)
    frameRate=1/180, # image acquisition rate in frames/sec
    max_lagtime=100, # intervals of frames out to which MSD is computed
         ):
    """
    Gets MSD from tracked cells
    """
    firstCell = None
    rcParams['figure.figsize']=5,3

    # track particles between frames
    t = tp.link(fb, maxDist, memory=maxMissFrame)
    print("Found %d "%t.shape[0])

    # keeps only trajectories that last for a given number of frames
    t1 = tp.filter_stubs(t,minFrame)
    print("Retained %d particles"%t1.shape[0])

    # correct for overall drifting motion
    d = tp.compute_drift(t1)
    tm = tp.subtract_drift(t1.copy(), d)
    print("Drift \n",d)

    if tm.shape[0] < 1:
      print("WARNING, nothing identified")
      return None

    # compute individual msd
    print("WARNING: commented out imsd")
    if False:
      im = tp.imsd(tm, pixelSize, frameRate, max_lagtime)
      plt.figure()
      plt.title('Individual MSD')
      plt.plot(im.index,im,color='tab:blue',alpha=0.1)
      if name is not None:
          plt.gcf().savefig('individualMSD_{}.png'.format(name) )
    else:
      im = False
    # compute ensemble msd
    em = tp.emsd(tm, pixelSize, frameRate, max_lagtime)

    # plot trajectories
    plt.figure()
    tp.plot_traj(t1)
    if name is not None:
        plt.gcf().savefig('traj_{}.png'.format(name) )

    # reformat MSD data
    time1 = []
    for i in em.index:
        time1.append(i)
    data = []
    for i in em.values:
        data.append(i)
    firstCell = np.array( data )

    print("Double check, data contains ALL particle trajectories? (%d)")



    plt.figure()
    plt.title('Ensemble MSD')
    plt.plot( time1,firstCell )
    plt.xlabel('sec')
    plt.ylabel('um**2')
    if name is not None:
        plt.gcf().savefig('ensembleMSD_{}.png'.format(name) )


    return time1,t1, data, im, em

def DoSave(time1,dataavg1,fileName="expt_noatp.csv"):
    """
    Saves input data as csv. Not sure if I need this
    """
    outData = np.zeros_like(np.outer(time1,[0,0]))
    outData[:,0]=time1
    outData[:,1]=dataavg1

    np.savetxt(fileName,outData)



def CalcDistances(t1):
    """
    Calculates distances covered by a particle over the entire simulation range
    t1 - traj object
    """
    particles = np.unique(t1['particle'])

    minFrames = 300
    dists = []
    xdists = []
    for particleNum in particles:
      subDF = t1.loc[t1['particle'] == particleNum]
      #subDF.loc[0]
      daMin = subDF['frame'].min()
      daMax = subDF['frame'].max()

      if daMin>0 or (daMax-daMin)<minFrames:
            next

      xi = subDF.loc[daMin].x
      yi = subDF.loc[daMin].y
      xf = subDF.loc[daMax].x
      yf = subDF.loc[daMax].y

      dist = (xf-xi)**2 + (yf-yi)**2
      dist = np.sqrt(dist)
      #dist = xf-xi
      dists.append( dist)

      xdists.append(xf - xi)

    return dists,xdists

def ProcessFrames(
  frames,
  downsampleRate = 0,
  crop = False,
  thresh=True
  ):
  """
  Downsample, crop and threshold image
  - frames: tiff file
  - downsampleRate: take every nth frame if >0
  - threshold: hardocded 
  """

  if downsampleRate>1:
    downsampled = frames[::downsampleRate,:,:]
  else:
    downsampled = frames


  if crop:
    cropped = downsampled[:,800:1100,800:1100]
    #plt.imshow(cropped[0,:,:])
  else:
    cropped = downsampled

  if thresh:
    threshed=np.zeros_like(cropped)
    threshed[np.where(cropped>105)]=255
  else:
    threshed=cropped


  print(np.shape(threshed))
  plt.imshow(threshed[-1,:,:],cmap="gray")
  plt.colorbar()
  return threshed

def doMSDFit(
  average_MSD_files,
  #fittingRange = [0,100]
  fittingRange=None
  ):
  """
  fits MSD to input csv files
  """
  
  for file, info in average_MSD_files.items():
      df = pd.read_csv(file)
      plt.errorbar(df['lagt'], df['msd'], yerr=np.std(df['msd']), label=f'{file}', color=info["color"])
      if fittingRange is not None:
        sub = df.iloc[0:fittingRange[1]]
      else:
        sub = df
      lagt = sub['lagt']
      msd = sub['msd']
      fit_coefficients = np.polyfit(lagt, msd, 1)
      linear_fit = np.poly1d(fit_coefficients)
      plt.plot(df['lagt'], linear_fit(df['lagt']), label=f'{file} - Linear Fit', linestyle='--', color=info["color"])


  plt.xlabel('lagt')
  plt.ylabel('msd')
  plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  plt.title('Average MSD with Standard Deviation and Linear Fit')
  print("slope/intercept",fit_coefficients)


  plt.savefig('averageMSD_all_121323.png', bbox_inches='tight')  # Adjust the file name and format as needed
  



