
"""
Script for validation routines

"""

def UnitTestGradient():
    """
    Test estimate of Javg for uncrowded domain. Flux should be the same everywhere in the domain, therefore Flux at some random point should equal that averaged over some region 
    """
    # each unit should have a gradient of -1 , no occlusions
    dims=[100,100]
    xs = 1-np.arange(dims[0])/dims[0]
    ys = np.ones(dims[1])
    data = np.outer(10*xs,ys) + np.reshape(0.01*np.random.rand(np.prod(dims)),dims)
    
    #data[30:60,30:60]=0
    
    gradc = np.diff(data,axis=0)
    #gradc[27:63,27:63]=0
    #
    D =1
    J = gradc * D
    J = J[5,2]
    print(J)
    
    areaFrac,Javg = CalcAverageFlux(data,D=1,xlims=[25,65])
    print(Javg)
    err = 1e-1
    plt.figure()
    plt.hist(gradc)
    #print('ss',Javg)
    assert np.abs(areaFrac-1)<err, "Unit test failed" # 1.0 --> no occlusions
    assert np.abs(Javg-J)<err, "Unit test failed" # 1.0 --> no occlusions
    print("PASS")

def UnitTestInclusion():
    """
    Test estimate of Javg for crowded domain (1 layer-like inclusion). Averaged flux should be proportional to the height-fraction occupied by the impenetrable layer
    """  
    # each unit should have a gradient of -1 , no occlusions
    dims=[100,100]
    xs = 1-np.arange(dims[0])/dims[0]
    ys = np.ones(dims[1])
    c0=10
    # adding noise 
    data = np.outer(c0*xs,ys) + 0.0*np.reshape(0.01*np.random.rand(np.prod(dims)),dims)
    
    w=30
    h=30
    data[30:(30+w),40:(40+h)]=0
    
    gradc = np.diff(data,axis=0)
    #gradc[27:63,27:63]=0
    #
    D =1
    J = gradc * D
    J = J[5,2]  # random pt far from inclusion that should have the correct linear conc. gradient 
    
    h0=dims[1]
    hFrac = (h0 - h)/h0 # analytic result for one square inclusion, based on media in parallel
    Janaly = J * hFrac
    print("Javg (Analytic)",Janaly)
    
    areaFrac,Javg = CalcAverageFlux(data,D=D,xlims=[30+1,30+w-1]) # avoiding margins by discarding 1 px on either side
    print("Javg (numerical)",Javg)
    err = 1e-1
    #plt.figure()
    #plt.hist(gradc)
    assert np.abs(Javg-Janaly)<err, "Unit test failed" # 1.0 --> no occlusions
    print("PASS")
    
UnitTestGradient()
UnitTestInclusion()




