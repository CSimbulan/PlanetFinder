import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import fits as f

fit = {}
biasFit = []
flatFit = []
targetFit = []
bias = None
flat = None
mask = None

def loadData(name, directory = ""): #Loads data
    return pf.open("{0}{1}".format(directory, name))

def getFileList(name, directory = ""):
    fileList = np.genfromtxt(directory + name, dtype='str', delimiter = "\t")
    fileInfo = []
    for row in fileList:
        fileInfo.append(filter(None, row.split('  ')))
    return np.array(fileInfo)


def setupFits(fileList):
    global bias
    global flat
    global mask
    loadBias()
    loadFlat()
    loadMask()
    setupbias = True
    setupflat = True
    if bias is not None:
        setupbias = False
    if flat is not None:
        setupflat = False
    
    for fi in fileList:
##        if (fi[2] == "FLAT" and not setupflat) or (fi[2] == "BIAS" and not setupbias):
##            continue
        tempf = getFit(fi)
        fit[fi[0][:-4]] = tempf
        if tempf.obstype == "FLAT":
            flatFit.append(tempf)
        elif tempf.obstype == "BIAS":
            biasFit.append(tempf)
        elif tempf.obstype == "TARGET":
            targetFit.append(tempf)

    if setupbias:
        bias = setupBias()
    if setupflat:
        if loadFlatUnNormalized() == None:
            setupFlatUnNormalized()
        mask = setupMask()
        flat = setupFlat()

    for t in targetFit:
        t.calibrate(mask, bias, flat)


def setupCentroids(fit = targetFit):
    last = None
    i = 0
    for t in fit:
        i += 1
        #print i
        t.quickCentroids(last)
        
        last = t.centroids


def calculateFit(x, y):
    '''
    Source from Shelly, Lab #2 - Astronomy Spectroscopy
    Calculates m and c for the fit
    '''
    ma = np.array([[np.sum(x**2), np.sum(x)],[np.sum(x), len(x)]])
    mc = np.array([[np.sum(x*y)], [np.sum(y)]])

    mai = np.linalg.inv(ma)
    md = np.dot(mai, mc)

    mfit = md[0,0]
    cfit = md[1, 0]
    return (mfit, cfit) #where m and c is from mx + c

def calculateFitError(x, y):
    # Return the error in both m and c
    m, c = calculateFit(x, y)
    N = len(x)
    std = (((N - 2.0)**-1.0)*sum((y - (m*x + c))**2))**0.5
    denominator = float(N*sum(x**2) -sum(x)**2)
    mstd = ((N*std**2)/denominator)**0.5
    cstd = ((sum(x**2)*std**2)/denominator)**0.5

    return (mstd, cstd)

def getNoiseGain(m, std):
    v = std**2.0
    k, vo = calculateFit(m, v)
    em, ec = calculateFitError(m, v)
    x = np.arange(min(m), max(m) + 1.0)
    y = vo + k*x
    
    plt.figure()
    plt.plot(m, v, '.', label="")
    plt.plot(x, y, 'r', label="Linear fit")
    plt.xlabel("Mean")
    plt.ylabel("Variance")
    plt.title("Mean and Variance relation")
    print "Noise = " + str(np.absolute(vo)**0.5) + "\nGain = " + str(k)
    print "Noise = {0} \nGain = {1}".format(ec**0.5, em)
    plt.legend(loc=4)
    plt.show()

def readoutNoise(fit = flatFit):
    global bias
    data = []
    l= (0, 0)
    for fl in flatFit:
        fl.loadData()
        l = fl.length()
        fl.data -= bias
        fl.data *= 30./fl.exptime

        data.append(fl.data)
    print "lala"
    median, std = getMedianStd(data)
    return (median, std)


                
        
def setupIntensity(fit = targetFit):

    import time

    start = time.time()

    circle = []
    annulus = []
    for i in range(30):
        circle.append(getCircle(i))
    circle = np.array(circle)

    outerCircle = getCircle(30)
    innerCircle = getCircle(25)

    x = []
    y = []
    for i in range(len(outerCircle[0])):
        add = True
        for j in range(len(innerCircle[0])):
            if outerCircle[0][i] == innerCircle[0][j] and outerCircle[1][i] == innerCircle[1][j]:
                add = False
                break
        if add:
            x.append(outerCircle[0][i])
            y.append(outerCircle[1][i])
    annulus = np.array([x, y])

    end = time.time()
    print end - start

    i = 0
    import time

    past = None
    start = time.time()
    for t in fit:
        t.quickIntensity(circle, annulus, past)
        past = t
    end = time.time()
    print end - start


def getCircle(r):
    theta = np.arange(0, 2*np.pi + 0.1, np.pi/(4*r + 0.1))
    xShell = np.around(r*np.cos(theta)).astype(int)
    yShell = np.around(r*np.sin(theta)).astype(int)
    x = np.array([])
    y = np.array([])
    for xi in range(np.min(xShell), np.max(xShell) + 1):
        i = np.where(xShell == xi)
        line = np.arange(np.min(yShell[i]), np.max(yShell[i]) + 1)
        x = np.concatenate((x, np.zeros(len(line)) + xi))
        y = np.concatenate((y, line))

    return np.array([x, y])
            
def loadBias():
    global bias
    bias = loadFit("bias.fit")
    return bias
    

def loadFlat():
    global flat
    flat = loadFit("flat.fit")
    return flat

def loadFlatUnNormalized():
    f = loadFit("flatUnNormalized.fit")
    return f

def loadMask():
    global mask
    mask = loadFit("mask.fit")
    return mask
    
    
    
    
def loadtxt(name, directory = ""):
    try: # textfile exists,
        f = np.loadtxt("{0}{1}".format(directory,name))
        return f
    except: # if not, create a new average text file
        return None

def loadFit(name):
    try: # textfile exists,
        return pf.getdata(f.DIR + name)
    except: # if not
        return None


def getFit(fileInfo):
    return f.fits(fileInfo)


def setupBias():
    data = []
    for b in biasFit:
        b.loadData()
        l = b.length()
        if l == (2501, 2148):
            b.data = b.data[130:-300]
        data.append(b.data)
    
    median = getMedian(data)
    savefit(median, "bias")
    return median


def getMedian(listData):
    i, j, k  = (len(listData[0]), len(listData[0][0]), len(listData))
    data = np.zeros((i,j))
    for row in range(i):
        for coloumn in range(j):
            l = []
            for f in range(k):
                l.append(listData[f][row][coloumn])
            data[row][coloumn] = np.median(l)
    return data

def getMedianStd(listData):
    i, j, k  = (len(listData[0]), len(listData[0][0]), len(listData))
    median = np.zeros((i,j))
    std = np.zeros((i,j))
    for row in range(i):
        for coloumn in range(j):
            l = []
            for f in range(k):
                l.append(listData[f][row][coloumn])
            median[row][coloumn] = np.median(l)
            std[row][coloumn] = np.std(l)
    return np.array([median, std])

def setupMask():
    flatData = loadFlatUnNormalized()
    l = (len(flatData), len(flatData[0]))
    mask = np.zeros(l)
    cutoff = np.mean(flatData) - np.std(flatData)
    for i in range(l[0]):
        for j in range(l[1]):
            if flatData[i][j] < cutoff:
                mask[i][j] = 0
            else:
                mask[i][j] = 1
    savefit(mask, "mask")
    return mask
        
    
def setupFlatUnNormalized():
    global bias
    data = []
    l= (0, 0)
    for fl in flatFit:
        fl.loadData()
        l = fl.length()
        """fl.data -= bias"""
        fl.data = np.subtract(fl.data, bias)
        fl.data *= 30./fl.exptime

        data.append(fl.data)
    flatMedian = getMedian(data)
    savefit(flatMedian, "flatUnNormalized")
    return flatMedian

def setupFlat():
    global mask
    global flat
    flatData = loadFlatUnNormalized()
    l = (len(mask), len(mask[0]))
    data = []
    for i in range(l[0]):
        for j in range(l[1]):
            if mask[i][j] == 1:
                data.append(flatData[i][j])
    median = np.median(data)
    flat = flatData/median
    
    savefit(flat, "flat")
    
    return flat


def saveCentroids(target = targetFit):
    i = 0
    for t in target:
        print i
        i += 1
        try:
            np.loadtxt(f.DIR + t.filename[:-4] + "_centroids.txt")
            continue
        except:
            savetxt(t.centroidsPos, "{0}_centroids".format(t.filename[:-4]))
        else:
            pass
        
    
def savetxt(data, name):
    np.savetxt("{0}{1}.txt".format(f.DIR, name), data)

def savefit(data, name):
    pf.writeto(f.DIR + name + ".fit", data, clobber=True )


def mask(data):
    mag = np.where(data > 1.2)
    

def plotLightCurve():
    I = []
    m = []
    for c in targetFit:
        I.append(c.target.intensity)
        m.append(c.I_r)
    I = np.array(I)
    m = np.array(m)
    r = I/m
    N = np.sum(r)/len(r)
    
    plt.plot(r/(N), '.')
    plt.show()

def lightCurve(targetNum):
    I = []
    I_e = []
    m = []
    t = []
    E = []
    
    refE = []
    refI = []
    out = 6
    if targetNum == 2 or targetNum == 3 or targetNum == 4 or targetNum == 1 or targetNum == 10:
        out = 5
    for rI in range(len(targetFit[0].centroids) - out):
        refI.append([])
        refE.append([])

    
    for c in targetFit:
        temp = []
        j = 0
        for i in range(len(c.centroids)):
            if i == targetNum:
                I.append(c.centroids[i].intensity)
                I_e.append(c.centroids[i].E)
            elif i != 2 and i != 3 and i != 4 and i!= 1 and i !=10:
                I_r = c.centroids[i].intensity
                E1 = c.centroids[i].E
                refI[j].append(I_r)
                refE[j].append(E1)
                j +=1
                
        
        t.append(c.time)
        
    temp = []

##    for i in range(len(refI)):
##        refI[i] = np.array(refI[i])/(np.sum(refI[i])/len(refI[i]))
##        refE[i] = np.array(refE[i])/(np.sum(refE[i])/len(refE[i]))

    for j in range(len(refI[0])):
        temp = []
        temp2 = []
        for i in range(len(refI)):
            temp.append(refI[i][j])
            temp2.append(refE[i][j])
        m.append(weightedAverage2(temp, temp2))
    
    I = np.array(I)
    m = np.array(m)
    t = np.array(t)
    I_e = np.array(I_e)
    
    N = np.median(I/m)

    plt.plot(t - t[0], I/m / N, '.b')
    tBinned = binningData(t, 8)[0]
    lBinned, eBinned = binningData(I/m / N, 8)
    plt.plot(tBinned - tBinned[0], lBinned, 'or')
    plt.errorbar(tBinned - tBinned[0], lBinned, yerr=eBinned, fmt='+r')

    plt.xlabel("Time (minutes)")
    plt.ylabel("Relative Flux")
    plt.title("Light Curve of Gj1214")

    ti1 = 75.4514
    ti2 = 90.0545
    te1 = 109.717
    te2 = 123.806

    tF =60*(te1 - ti2)
    tT = 60*(te2 - ti1)
    #print tT
    #print tF
    #tT = 52.753*60
    #tF = 28.7*60
    i1 = np.where(t <=75.4514 + t[0])[0]
    i2 = np.where(t >= 123.806 + t[0])[0]
    i = np.concatenate((i1,i2))
    F = np.median(((I/m)/N)[i])
    Fstd = np.std(((I/m)/N)[i])


    plt.plot(t[i1] - t[0], np.zeros(len(t[i1])) + F, '-k', linewidth=2)
    plt.plot(t[i2] - t[0], np.zeros(len(t[i2])) + F, '-k', linewidth=2 )


    i1 = np.where(t >= ti2 + t[0])[0]
    i2 = np.where(t[i1] <= te1 + t[0])[0] + np.min(i1)
    dF = np.mean(((I/m)/N)[i2])

    tLine=np.arange(dF, F, 0.0001 )
    plt.plot(ti1 + np.zeros(len(tLine)), tLine, '--g', linewidth=2)
    plt.plot(te2 + np.zeros(len(tLine)), tLine, '--g', linewidth=2)

    plt.plot(ti2 + np.zeros(len(tLine)), tLine, '--m', linewidth=2)
    plt.plot(te1 + np.zeros(len(tLine)), tLine, '--m', linewidth=2)

    plt.plot(t[i2] - t[0], np.zeros(len(t[i2])) + dF, '-k', linewidth=2)
    plt.show()
    
    dFstd = np.std(((I/m)/N)[i2])
    
    delta = 1 - dF/F #DELTA
    deltaE = delta*((dFstd/dF)**2 + (Fstd/F)**2)**0.5 #DELTA ERROR
    print "delta = dF/F = {0} +/- {1}".format(delta, deltaE)

    import uncertainties as u
    from uncertainties import unumpy
    delta = u.ufloat(delta, deltaE)
    print "delta = dF/F = {0}".format(delta)

    
    Rsun = 695800000
    Rstar = u.ufloat(0.216*(Rsun),0.012*Rsun )

    RJupiter = (69911000)
    Rexo = (Rstar*(delta)**0.5)/RJupiter

    print "Radius = {0}".format(Rexo)


    #Rp = u.ufloat(Rexo, RexoE)*RJupiter
    
    G = 6.67384e-11
    Msun = 1.989e30
    Ms = u.ufloat(0.15*Msun, 0.011*Msun)
    MJupiter = 1.898e27
    Mexo = u.ufloat(0.0203253311519,0.00311487273071)*MJupiter
    rho = Mexo/(4./3 * np.pi * (Rexo*RJupiter)**3)
    print "Density = {0} g/cm^3".format(1000*rho/100**3)
    P = Ms*G*np.pi*(tT**2 - tF**2)**(1.5) / (32 * Rstar**3 * delta**(0.75))
    print "Period = {0} Days".format(P/ 60 / 60 / 24)

    g = G*Mexo/(Rexo*RJupiter)**2
    print "g = {0} cm/s^2".format(g*100)

    
    a = (P**2 * G * Ms / (4 * np.pi**2))**(1./3)
    b = ((    (1 - delta**0.5)**2    -   (tF/tT)**2 * (1 + delta**0.5)**2 )/(1 - (tF/tT)**2))**0.5

    au = 149597870700
    print "a = {0} AU".format(a/au)
    print "b = {0} ".format(b)
    #b = 0.71
    
    i = u.unumpy.arccos(b*Rstar/a)
    print "Inclination = {0} degrees".format(i*180/np.pi)
    T = P/np.pi * u.unumpy.arcsin(((Rstar + Rexo*RJupiter)**2 - b**2)**.5 / a)
    print "Transit Time = {0} minutes".format(T/60.)
    print
    print

    p = Rexo * RJupiter / Rstar
    print p
    t1 = np.arange(-T.n/2., T.n/2. + 1, 1)

    
    i = u.ufloat(1,0)*i
    
    z = (a.n / Rstar.n) * ((np.sin(t1*2*np.pi/P.n))**2 + (np.cos(i.n)*np.cos(t1*2*np.pi/P.n))**2)**.5
    from occultquad import *
    x=occultquad(z, 1.0058, -0.1264, p.n)

    plt.xlabel("Relative Time (s)")
    plt.ylabel("Relative Flux")
    plt.title("limb darkening")
    plt.plot(t1, x[0], label="With limb darkening")
    plt.plot(t1, x[1], label ="Without limb darkening")
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.plot(t - t[0], I/m / N, '.b')
    plt.plot(tBinned - tBinned[0], lBinned, 'or', label="Observed Curve")
    plt.errorbar(tBinned - tBinned[0], lBinned, yerr=eBinned, fmt='+r')
    
    plt.plot(t1/60 + 99., x[0], 'm', label="With Limb darkening")
    plt.plot(t1/60 + 99., x[1], 'g',label="Without limb darkening")
    plt.xlabel("Time (minutes)")
    plt.ylabel("Relative Flux")
    plt.title("Relative transit flux")
    plt.legend(loc="upper center")
    plt.show()

    return x

    
def weightedAverage(values):
    values = np.array(values)
    return np.sum(values**-1)/np.sum(values**-2.)

def averageError(values):
    values = np.array(values)
    return np.std(values)

def weightedAverage2(flux, error):
    flux = np.array(flux)
    error = np.array(error)
    w = error**(-2)
    return np.sum(w*flux)/np.sum(w)

def binningData(data, l):
    v = []
    e = []
    t = []
    for i in range(len(data)):
        t.append(data[i])
        if len(t) == l:
            v.append(weightedAverage(t))
            e.append(averageError(t))
            t = []
    return np.array([v, e])


def saveIntensity(fit = targetFit):
    for f in fit:
        pass
        
    
    
def quickLightCurve():
    i = len(targetFit[0].centroids)
    
    for j in range(0,i):
        plt.subplot(3, 4, j)
        print j
        lightCurve(j)
    plt.show()
    #plt.legend()
    #plt.show()
    
        
	
    

if __name__ == "__main__":
    print "Getting filelist"
    fileList = getFileList("filelist", f.DIR)
    print "Setting up fits"
    setupFits(fileList)
    print "Setting up centroids"
    setupCentroids()
    print "Calculating intensity"
    setupIntensity()
    #saveCentroids()
