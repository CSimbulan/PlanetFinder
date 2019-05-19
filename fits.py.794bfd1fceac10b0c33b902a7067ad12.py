import numpy as np
import astropy.io.fits as pf
import urllib as url
import matplotlib.pyplot as plt
from star import *

filename, targetname, obstype, filt, exptime = None, None, None, None, None
DIR = "./ACAM_DATA/"


READOUTNOISE = 33.0968452445
READOUTERROR = 1.43800861993 


        
class fits:
    header = None
    data = None

    peaks = None
    centroidsPos = None
    centroids = []

    target = None

    I_r = 0
    time = 0

    
    
    def __init__(self, info):
        self.filename, self.targetname, self.obstype, self.filt, self.exptime = info
        self.exptime = float(self.exptime)
        self.defineFit(self.filename)
        self.time = self.setTime()
        
    def defineFit(self, filename):
        fits1 = DIR + filename
        self.header = pf.getheader(fits1)

    def setTime(self):
        time = self.header["ST"][1:]
        h,mi,s = int(time[0:2]), int(time[3:5]), float(time[6:])
        return ((h)*60 + mi) + s/60.
        
    def length(self):
        """
        if self.data == None:
            return None
        else:
        """
        return (len(self.data), len(self.data[0]))

    def loadData(self):
        if self.data == None:
            fits1 = DIR + self.filename
            self.data = pf.getdata(fits1)

    def normalize(self, bias): #Normalize the data with respect to its mean
        i, j = self.length()
        if self.obstype == "FLAT" and (("HISTORY" not in self.header) or ("HISTORY" in self.header and "NORMALIZED" not in self.header["HISTORY"])):
            self.data /= np.mean(self.data - bias[(i, j)])
            #self.saveFit(msg = "NORMALIZED")

    
    def calibrate(self, mask, bias, flat):
        if self.obstype != "TARGET":
            return
        #If statements checks whether its already calibrated or not
        if ("HISTORY" not in self.header) or ("HISTORY" in self.header and "CALIBRATED" not in self.header["HISTORY"]):
            self.loadData()
            l = self.length()
            self.data = (self.data - bias)/flat
            self.data = np.nan_to_num(self.data)
            self.data *= mask ## Remove this if you dont want to save it masked
            self.saveFit(msg = "CALIBRATED")

    def saveFit(self, msg="updated"):
        self.header["HISTORY"] = msg
        pf.writeto(DIR + self.filename, self.data, header = self.header,clobber=True )

    def quickCentroids(self, lastCentroids = None, cut = 5000, radius = 3):
        self.loadData()
        self.centroids = []
        test = False
        try:
            self.centroidsPos = np.loadtxt(DIR + self.filename[:-4] + "_centroids.txt")
        except:
            if lastCentroids != None:
                test = True
                self.findStars(lastCentroids, 5)
            else:
                
                self.cutOff(cut)
                self.findCentroids(radius)
        else:
            pass
        for i in range(len(self.centroidsPos[0])):
            p = (self.centroidsPos[0][i], self.centroidsPos[1][i])
            s = star(p)
            self.centroids.append(s)
            if s.target:
                self.target = s
        self.centroids = np.array(self.centroids)
            

    def findStars(self, ref, r = APERTURE_RADIUS):
        c = []
        self.centroidsPos = [[],[]]
        for re in ref:
            p = np.around(np.array(re.location)).astype(int)
            maxI = 0
            maxPos = p
            for x in range(p[0] - r, p[0] + r + 1):
                for y in range(p[1] - r, p[1] + r + 1):
                    if self.data[y][x] > maxI:
                        maxI = self.data[y][x]
                        maxPos = (x, y)

            ce = self.determineCentroids(maxPos, 5)
            c.append(ce)
            self.centroidsPos[0].append(ce[0])
            self.centroidsPos[1].append(ce[1])
        
        #self.centroidsPos = np.array(c)
            
            
                
        

    def plot(self):
        plt.imshow(self.data)
        theta = np.arange(0, 2*np.pi + 0.1, 0.1)
        for c in self.centroids:
            ri = c.aperture_radius
            ris = APERTURE_INNERSKY
            ros = APERTURE_OUTERSKY
            plt.plot(ri*np.cos(theta) + c.location[0], ri*np.sin(theta) + c.location[1], 'g')
            plt.plot(ris*np.cos(theta) + c.location[0], ris*np.sin(theta) + c.location[1], 'g')
            plt.plot(ros*np.cos(theta) + c.location[0], ros*np.sin(theta) + c.location[1], 'g')
            #plt.plot(c.location[0], c.location[1], '.r')

        plt.xlim(0, len(self.data[0]))
        plt.ylim(0, len(self.data))
        plt.xlabel("x [Pixels]")
        plt.ylabel("y [Pixels]")
        plt.title("{0}".format(self.filename[:-4]))
        plt.show()           

    def cutOff(self,s):
        data = self.data
        bright = np.where(data > s)
        self.peaks = bright[::-1]
        return bright



##    def getIntensity(self, pos, ri, ro):
##        ri = int(ri)
##        ro = int(ro)
##        pos = np.around(np.array(pos)).astype(int)
##        intensity = []
##        for rx in range(0, ro + 1):
##            for ry in range(0, ro - np.abs(rx) + 1):
##                if rx < ri and ry < ri:
##                    continue
##                elif ry == 0 and rx == 0:
##                    intensity.append(self.data[pos[1]][pos[0]])
##                elif ry == 0 and rx>=ri:
##                    intensity.append(self.data[pos[1] + rx][pos[0]])
##                    intensity.append(self.data[pos[1] - rx][pos[0]])
##                elif rx == 0 and ry >= ri:
##                    intensity.append(self.data[pos[1]][pos[0] + ry])
##                    intensity.append(self.data[pos[1]][pos[0] - ry])
##                else:
##                    intensity.append(self.data[pos[1] + rx][pos[0] + ry])
##                    intensity.append(self.data[pos[1] - rx][pos[0] + ry])
##                    intensity.append(self.data[pos[1] + rx][pos[0] - ry])
##                    intensity.append(self.data[pos[1] - rx][pos[0] - ry])
##                
##                    
##        return intensity

##    def getIntensity(self, pos, ri, ro = -1):
##        ri = int(ri)
##        ro = int(ro)
##        pos = np.around(np.array(pos)).astype(int)[::-1]
##        fullCircle = getCircle(pos, ri)
##        removeCircle = np.array([[],[]])
##        if ro != -1:
##            removeCircle = getCircle(pos, ro)
##
##        I = []
##        for i in range(len(fullCircle[0])):
##            add = True
##            for j in range(len(removeCircle[0])):
##                if fullCircle[0][i] == removeCircle[0][j] and fullCircle[1][i] == removeCircle[1][j]:
##                    add = False
##                    break
##            if add:
##                I.append(self.data[fullCircle[0][i]][fullCircle[1][i]])
##        return np.array(I)

    def getIntensity(self, circle):
        intensity = []
        x, y = circle

        for i in range(len(x)):
            intensity.append(self.data[x[i]][y[i]])

        return np.array(intensity)


    def getRadius(self, ri, ro):
        ri = int(ri)
        ro = int(ro)
        r = []
        for rx in range(0, ro + 1):
            for ry in range(0, ro - np.abs(rx) + 1):
                if rx < ri and ry < ri:
                    continue
                elif ry == 0 and rx == 0:
                    r.append(0)
                elif ry == 0 and rx>=ri:
                    r.append(D(rx, 0))
                    r.append(D(rx, 0))
                elif rx == 0 and ry >= ri:
                    r.append(D(0, ry))
                    r.append(D(0, ry))
                else:
                    r.append(D(rx, ry))
                    r.append(D(-rx, ry))
                    r.append(D(rx, -ry))
                    r.append(D(-rx, -ry))

        return np.array(r)

    def std(self, x):
        x = np.array(x)
        #return (np.sum((x - np.median(x))**2)/len(x))**0.5
        #return np.mean((x - np.median(x))**2.)**0.5
        medianX = np.median(x)
        medianX2 = np.median(x**2)
        return (medianX2 - medianX**2)**0.5

                    
    def quickIntensity(self, aperture, annulus, past = None):
        if past == None:
            for s in self.centroids:
                Fn = []
                pos = np.around(s.location).astype(int)
                ann = (annulus[0] + pos[1], annulus[1] + pos[0])
                bgNoise = self.getIntensity(ann)
                bgMedian = np.median(bgNoise)
                bgstd = self.std(bgNoise)
                E = []
                for i in range(0,len(aperture)):
                    ap = (aperture[i][0] + pos[1], aperture[i][1] + pos[0])
                    tempi = self.getIntensity(ap)
                    Fn.append(np.sum(tempi - bgMedian))
                    E.append(((len(tempi)*(1 + len(bgNoise)**-1.))*(bgstd**2) + Fn[-1])**0.5)
                Fn = np.array(Fn)
                E = np.array(E)


                s.growthCurve = Fn

                #plt.plot(Fn/E, '--')
                #plt.show()
                

                

                #E = (Fn + N1*(expB + vRO) + N1*(expB + vRO)/N23)**0.5
                s.E = E
                snr = Fn/E
                s.snr = snr
                m = np.max(snr)
                def find90Max(v):
                    max90 = 0.975*(np.max(v) - np.min(v)) + np.min(v)
                    shortest = -1
                    lim = np.where(v == np.max(4))[0]
                    for s in snr:
                        if s> max90 :
                            break
                        if shortest == -1 or shortest > np.abs(max90 - s ):
                            shortest = s
                    return np.where(shortest == v)
                r_aperture = find90Max(snr)[0]
                #print r_aperture
                #r_aperture = np.where(snr == m)
                s.aperture_radius = r_aperture[0]
                
                s.intensity = Fn[r_aperture][0]
                s.E = E[r_aperture][0]

                #plt.figure()
                #plt.plot(Fn/E, '-k')
    ##                plt.xlabel("Aperure Radius")
    ##                plt.ylabel("Signal to Noise Ratio")
    ##                plt.title("Signal to Noise Ratio")
               # plt.show()

        else:
            for i in range(len(self.centroids)):
                s = self.centroids[i]
                s.aperture_radius = past.centroids[i].aperture_radius


                Fn = []
                pos = np.around(s.location).astype(int)
                ann = (annulus[0] + pos[1], annulus[1] + pos[0])
                ann = np.array(ann).astype(int)
                bgNoise = self.getIntensity(ann)
                bgMedian = np.median(bgNoise)
                bgstd = self.std(bgNoise)
                E = []
                area = (aperture[s.aperture_radius][0] + pos[1], aperture[s.aperture_radius][1] + pos[0])
                tempi = self.getIntensity(area)
                Fn = np.sum(tempi - bgMedian)
                E = ((len(tempi))*(bgstd**2) + Fn)**0.5
              
                Fn = np.array(Fn)
                E = np.array(E)

                s.E = E
                s.intensity = np.sum(Fn)

    def weightedAverage(self, flux, error):
        values = np.array(values)
        w = error**(-2)
        return np.sum(w*flux)/np.sum(w)
                
    def findCentroids(self,radius = 3):
        peaks = self.peaks[::-1]
        data = self.data
        coord = toCoordinates(peaks) # Get in (x, y) coordinates
        groups = grouping(coord) ##Get the groupings

        eX = []
        eY = []
        for g in groups:
            if len(g) == 0: continue 
            g = addToGroup(g, radius)
            I = 0
            Ii = []
            for p in g: #gets intensity for each coordinate
                x, y = p
                Ii.append(data[y][x])

            I = sum(Ii)
            Ix = {}
            Iy = {}
            for i in range(len(g)):
                if g[i][0] not in Ix.keys():
                    Ix[g[i][0]] = Ii[i]
                else:
                    Ix[g[i][0]] += Ii[i]
                if g[i][1] not in Iy.keys():
                    Iy[g[i][1]] = Ii[i]
                else:
                    Iy[g[i][1]] += Ii[i]
            ex = 0
            ey = 0
            for ix in Ix.keys():
                ex += ix*Ix[ix]
            for iy in Iy.keys():
                ey += iy*Iy[iy]
            eX.append(ex/I)
            eY.append(ey/I)
            
        self.centroidsPos = np.array([eX, eY])
        return np.array([eX, eY])

    def determineCentroids(self, pos, radius):
        g = addToGroup(np.array([pos]), radius)
        I = 0
        Ii = []
        for p in g: #gets intensity for each coordinate
                x, y = p
                Ii.append(self.data[y][x])
        I = sum(Ii)
        Ix = {}
        Iy = {}
        for i in range(len(g)):
            if g[i][0] not in Ix.keys():
                Ix[g[i][0]] = Ii[i]
            else:
                Ix[g[i][0]] += Ii[i]
            if g[i][1] not in Iy.keys():
                Iy[g[i][1]] = Ii[i]
            else:
                Iy[g[i][1]] += Ii[i]
        ex = 0
        ey = 0
        for ix in Ix.keys():
            ex += ix*Ix[ix]
        for iy in Iy.keys():
            ey += iy*Iy[iy]
        eX = (ex/I)
        eY = (ey/I)
        return np.array([eX, eY])
        

    def loadMask(self):
        global mask
        mask = loadFit("mask.fit")
        return mask
    
    def loadFit(self, name):
        try: # textfile exists,
            return pf.getdata(DIR + name)
        except: # if not
            return None
        

    def toDegrees(self):
        ras = self.header['ra'][1:]
        des = self.header['dec']
        radeg = 15*(float(ras[0:2]) + float(ras[3:5])/60. + float(ras[6:])/3600.)
        dsgn = np.sign(float(des[0:3]))
        dedeg = float(des[0:3]) + dsgn*float(des[4:6])/60. + dsgn*float(des[7:])/3600.
        fovam = 10. # size of square search field in arc min
        return (radeg, dedeg, fovam)



    
def removearray(L,arr):
    '''
    http://stackoverflow.com/questions/3157374/how-do-you-remove-a-numpy-array-from-a-list-of-numpy-arrays
    '''
    ind = 0
    size = len(L)
    while ind != size and not np.array_equal(L[ind],arr):
        ind += 1
    if ind != size:
        L.pop(ind)
    else:
	#raise ValueError('array not found in list.')
        pass
    return L        


def toCoordinates(points):
    points = np.array([points[1], points[0]])
    return np.array(tuple(zip(*points)))

def addToGroup(group, radius):
    newGroup = np.array(group)
    if radius < 0: return np.array(newGroup)
    for p in group:
        newGroup = addAdj(p, newGroup)
    newGroup = addToGroup(newGroup, radius - 1)
    return np.array(newGroup)

def addAdj(p, group):
    newGroup = list(group)
    adj = np.array([[p[0] + 1, p[1]], [p[0] - 1, p[1]], [p[0], p[1] + 1],[p[0], p[1] - 1]])
    for a in adj:
        if a[0] >= 2071 or a[1] >= 2148 or a[0] < 0 or a[1] < 0: continue
        b = False
        for i in newGroup:
            if (a == i).all():
                b = True
                break
        if b == False:
            newGroup.append(a)
    return np.array(newGroup)
                

def grouping(data):
    unGrouped = list(data)
    groups = []
    for d in data:
        #print (d == unGrouped).any(0).all()
        if len(unGrouped) == 0 or not(d == unGrouped).any(0).all():
            continue
        g = group(d, data)
        for i in g:
            removearray(unGrouped, i)
        groups.append(g)
        
    #for uG in unGrouped:
    #    groups.append([uG])
    return groups
    
def group(point, data):
    nearData = list(data)
    for d in data:
        if pixelDistance(point, d) >20:
            removearray(nearData, d)
    adjPoints = []
    star = [point]
    for d in nearData:
        if areAdjacent(point, d):
            adjPoints.append(d)
    star = adjPoints
    for a in adjPoints:
        removearray(nearData, a)
    for a in adjPoints:
        gr = group(a, nearData)
        for g in gr:
            removearray(nearData, g)
            star.append(g)
    return star

def areAdjacent(point1, point2):
    if (point2 == (point1[0] + 1, point1[1])).all(): return True
    if (point2 == (point1[0] - 1, point1[1])).all(): return True
    if (point2 == (point1[0], point1[1] + 1)).all(): return True
    if (point2 == (point1[0], point1[1] - 1)).all(): return True
    return False

def pixelDistance(P1, P2):
    px = P2[0] - P1[0]
    py = P2[1] - P1[1]
    return (px**2 + py**2)**0.5


def D(x,y):
    return (x**2 + y**2)**2

def getPosTest(pos, ri, ro):
        ri = int(ri)
        ro = int(ro)
        pos = np.around(np.array(pos)).astype(int)
        x = []
        y = []
        for rx in range(0, ro + 1):
            for ry in range(0, ro - np.abs(rx) + 1):
                if rx < ri and ry < ri:
                    continue
                elif ry == 0 and rx == 0:
		    x.append(pos[1])
		    y.append(pos[0])
                elif ry == 0 and rx>=ri:
		    x.append(pos[1] + rx)
		    y.append(pos[0])
		    x.append(pos[1] - rx)
		    y.append(pos[0])
                elif rx == 0 and ry >= ri:
		    x.append(pos[1])
		    y.append(pos[0] + ry)
		    x.append(pos[1])
		    y.append(pos[0] - ry)

                else:
		    x.append(pos[1] + rx)
		    y.append(pos[0] + ry)
		    x.append(pos[1] + rx)
		    y.append(pos[0] - ry)
		    x.append(pos[1] - rx)
		    y.append(pos[0] + ry)
		    x.append(pos[1] - rx)
		    y.append(pos[0] - ry)
	return np.array([x, y])

def getCircle(pos, r):
    theta = np.arange(0, 2*np.pi + 0.1, np.pi/(3*r + 0.1))
    xShell = np.around(r*np.cos(theta) + pos[0]).astype(int)
    yShell = np.around(r*np.sin(theta) + pos[1]).astype(int)
    x = np.array([])
    y = np.array([])
    for xi in range(np.min(xShell), np.max(xShell) + 1):
        i = np.where(xShell == xi)
        line = np.arange(np.min(yShell[i]), np.max(yShell[i]) + 1)
        x = np.concatenate((x, np.zeros(len(line)) + xi))
        y = np.concatenate((y, line))

    return np.array([x, y])

	
    


