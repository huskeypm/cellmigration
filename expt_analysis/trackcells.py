import matplotlib.pylab as plt
import numpy as np
import trackpy as tp 

# get tracking params from a certain frame
def TrialParams(frames,
              refFrame=0, # reference frame to find paticles
              diameter=41, # diameter of the particle in pixel, must be an odd number
              minmass=1e2, # total brightness of a particle: a threshold to discard floating particles
             ):

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

# track particles using params determined from getParams()
def DoTracking(frames,diameter=41,minmass=1e2):
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



