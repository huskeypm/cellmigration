"""
This program allows one to evolve the states of particles 
0 -> 1 -> 2
^         |
+---------+
using a gillespie algorithm 

"""
import numpy as np
import matplotlib.pylab as plt

# state, t0, k01, k12
IDXSTATE=0
IDXT0=1
IDXK01=2
IDXK12=3
IDXK20=4
IDXCA =5 # contact with A (0,1)
IDXCB =6  # contact with A (0,1)
# state - current state (0,1,2)
# t0 - last time state flip
# kij - rate to flip from i to j


class CellSystem():
  def __init__(self,nCells=10):
    self.stateMatrix = np.zeros([nCells,7])
    self.stateMatrix[:,IDXSTATE]=0
    self.stateMatrix[:,IDXT0   ]=0 # for now np.array(-100*np.random.randn(nCells),int) # start before 0
    self.stateMatrix[:,IDXK01  ]=1e-2
    self.stateMatrix[:,IDXK12  ]=1e-3
    self.stateMatrix[:,IDXK20  ]=1e-4
    self.stateMatrix[:,IDXCA   ]=0.    
    self.stateMatrix[:,IDXCB   ]=0.     
    self.nCells = nCells

  def Update(self,values,idx=IDXSTATE):
    """
    updates statematrix with current set of values    
    """
    self.stateMatrix[:,  idx   ]=values

  def UpdateK01(self,values):
    self.Update(values,idx=IDXK01)

  def UpdateK12(self,values):
    self.Update(values,idx=IDXK12)

  def UpdateK20(self,values):
    self.Update(values,idx=IDXK20)

  def UpdateContactsA(self,values):
    self.Update(values,idx=IDXCA)

  def UpdateContactsB(self,values):
    self.Update(values,idx=IDXCB)

  def PrintStates(self):
    print(self.stateMatrix[:,IDXSTATE])

  def SumStates(self):
    states = {'0':0,'1':0,'2':0}
    for key in states.keys():
      i = int(key)
      #states['%d'%i] np.count_nonzero(states[i,:] == i)
      states[key] = np.count_nonzero(self.stateMatrix[:,IDXSTATE] == i)
    return states

  def EvalTransitions(self,t):
        """
        uses gillespie algorithm to advance states
        """ 
        nCells = self.nCells
        stateMatrix = self.stateMatrix
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
        # is in state 0, Poccur is True, and in contact with CA
        idx0=np.where( (stateMatrix[:,IDXSTATE]==0) & (Poccur01) & (stateMatrix[:,IDXCA]==1))   
        # is in state 1, Poccur is True, and in contact with CB
        idx1=np.where( (stateMatrix[:,IDXSTATE]==1) & (Poccur12) & (stateMatrix[:,IDXCB]==1))
        # is in state 2, Poccur is True, 
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


def Iterator(cs): # ,stateMatrix,nCells):
  """
  cs - cell system 
  """
  stateMatrix = cs.stateMatrix
  nCells = cs.nCells
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
  
      cs.EvalTransitions(t) # ,nCells,stateMatrix)
  
      # update 
      states[i,:]=stateMatrix[:,IDXSTATE]
      #print(states[i,:])
      
      # inefficient
      stateSums = cs.SumStates()
      isZero[i] = stateSums['0'] #(states[i,:] == 0)
      isOne[i]  = stateSums['1'] #p.count_nonzero(states[i,:] == 1)
      isTwo[i]  = stateSums['2'] #p.count_nonzero(states[i,:] == 2)
      
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

#!/usr/bin/env python
import sys
##################################
#
# Revisions
#       10.08.10 inception
#
##################################

#
# ROUTINE  
#
def doit(fileIn):
  nCells = 5
  cs = CellSystem(nCells=nCells)

  cs.UpdateContactsA( np.ones(nCells) ) 
  cs.UpdateContactsB( np.ones(nCells) ) 

  t =0 
  cs.EvalTransitions(t) #,nCells,cs.stateMatrix) 

  Iterator(cs)


#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -validation" % (scriptName)
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
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  #fileIn= sys.argv[1]
  #if(len(sys.argv)==3):
  #  1
  #  #print "arg"

  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-validation"):
      doit("SDF")      
      quit()
  





  raise RuntimeError("Arguments not understood")




