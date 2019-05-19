import numpy as np
import matplotlib.pyplot as plt


TARGETLOC = (1417, 826) #target likely in this position
TARGETRADIUS = 50 #Check in this radius
APERTURE_RADIUS = 5.
APERTURE_INNERSKY = APERTURE_RADIUS*5
APERTURE_OUTERSKY = APERTURE_RADIUS*6

        
class star:
    location = (None, None)
    target = False

    growthCurve = None
    snr = None

    aperture_radius = APERTURE_RADIUS
    
    def __init__(self, location):
       self.location = location
       self.isTarget()


    def isTarget(self):
        d = pixelDistance(self.location, TARGETLOC)
        if d < TARGETRADIUS:
            self.target = True


    def setIntensity(self):
        self.Io_median = np.median(self.intensity_outersky)
        self.Io_std = np.std(self.intensity_outersky)
        self.intensity_radius -= self.Io_median

        self.I_r = np.sum(self.intensity_radius)

        
def pixelDistance(P1, P2):
    px = P2[0] - P1[0]
    py = P2[1] - P1[1]
    return (px**2 + py**2)**0.5

