import numpy as np 

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
        

