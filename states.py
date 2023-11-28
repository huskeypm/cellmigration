"""
This program allows one to evolve the states of particles 
0 -> 1 -> 2
^         |
+---------+
using a gillespie algorithm 

"""
import numpy as np
import matplotlib.pylab as plt

nCells = 20

# state, t0, k01, k12
IDXSTATE=0
IDXT0=1
IDXK01=2
IDXK12=3
IDXK20=4
# state - current state (0,1,2)
# t0 - last time state flip
# kij - rate to flip from i to j
stateMatrix = np.zeros([nCells,5])
stateMatrix[:,IDXSTATE]=0
stateMatrix[:,IDXT0   ]=0 # for now np.array(-100*np.random.randn(nCells),int) # start before 0
stateMatrix[:,IDXK01  ]=1e-2
stateMatrix[:,IDXK12  ]=1e-3
stateMatrix[:,IDXK20  ]=1e-4




def Iterator(stateMatrix):
  # reset 
  #stateMatrix[:,IDXT0   ]=0
  #stateMatrix[:,IDXSTATE]=0
  #np.random.seed(1)
  #r=np.random.rand(nCells)
  
  steps = 50
  interval = 100
  states=np.zeros([steps,nCells])
  isZero= np.zeros(steps)
  isOne = np.zeros(steps)
  isTwo = np.zeros(steps)
  times = np.arange(steps)*interval
  for i in range(steps):
      t = i*interval
      # Prob not occur
      #print(t-stateMatrix[:,IDXT0   ])
      def EvalTransitions(t,stateMatrix):
        w01 = np.exp(-(t-stateMatrix[:,IDXT0   ])*stateMatrix[:,IDXK01  ])
        w12 = np.exp(-(t-stateMatrix[:,IDXT0   ])*stateMatrix[:,IDXK12  ])  
        w20 = np.exp(-(t-stateMatrix[:,IDXT0   ])*stateMatrix[:,IDXK20  ])  
        # Prob will occur  (1-not)
        Poccur01 = np.array(np.random.rand(nCells) < 1-w01,int)
        Poccur12 = np.array(np.random.rand(nCells) < 1-w12,int)
        Poccur20 = np.array(np.random.rand(nCells) < 1-w20,int)
        
        
        # TESTING 
        #print(t,'r',r)
        #Poccur01 = r < (1-w01)
        #Poccur12 = r < (1-w12)
        #Poccur20 = r < (1-w20)
    
        # find which changed 
        idx0=np.where( (stateMatrix[:,IDXSTATE]==0) & (Poccur01))   
        idx1=np.where( (stateMatrix[:,IDXSTATE]==1) & (Poccur12))
        idx2=np.where( (stateMatrix[:,IDXSTATE]==2) & (Poccur20))
        # 0->1    
        stateMatrix[idx0,IDXSTATE]= 1
        # 1->2 
        stateMatrix[idx1,IDXSTATE]= 2
        # 2->0
        stateMatrix[idx2,IDXSTATE]= 0
        
        #print(idx0[0],idx1[0],idx2[0])
        #nChanged = np.shape(idx0[0])[0] +np.shape(idx1[0])[0] + np.shape(idx2[0])[0]
        #print(nChanged)
        stateMatrix[idx0,IDXT0   ]=t
        stateMatrix[idx1,IDXT0   ]=t
        stateMatrix[idx2,IDXT0   ]=t
  
      EvalTransitions(t,stateMatrix)
  
      # update 
      states[i,:]=stateMatrix[:,IDXSTATE]
      #print(states[i,:])
      
      # inefficient
      isZero[i] = np.count_nonzero(states[i,:] == 0)
      isOne[i]  = np.count_nonzero(states[i,:] == 1)
      isTwo[i]  = np.count_nonzero(states[i,:] == 2)
      
  #print(isZero)   
  #print(isOne)   
  #print(isTwo)   
  plt.figure()
  plt.plot(times,isZero, label=0)
  plt.plot(times,isOne,label=1)
  plt.plot(times,isTwo,label=2)
#plt.plot(times,isZero+isOne+isTwo)
  plt.legend(loc=0)
  plt.gcf().savefig('test.png')

Iterator(stateMatrix)
