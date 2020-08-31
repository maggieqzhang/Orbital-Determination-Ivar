# OD Photometry

import math
from astropy.io import fits
import matplotlib.pyplot as plt

def fractionalPhotometry(image):
    im = fits.getdata(image)
    
    print("Fractional Photometry")
    astX = int(input("Enter the x coordinate of the asteroid"))
    astY = int(input("Enter the y coordinate of the asteroid"))
    astRadius = int(input("Enter the radius of the circle"))
    asteroid = im[astX-astRadius:astX+astRadius+1, astY-astRadius:astY+astRadius+1]
    print(asteroid)
     
    starSkyCount = asteroid.sum()
    skyInner = int(input("Enter inner annulus radius"))
    skyOuter = int(input("Enter outer annulus radius"))
    sky = calcSky(astX, astY, skyInner, skyOuter)
    
    signal = starSkyCount-avgSky
    magnitude = -2.5*math.log10(signal)
    constant = mag-magnitude
    
    astHalf = int(astLength/2)
    astObj = im[astY-astHalf-2:astY+astHalf+3, astX-astHalf-2:astX+astHalf+3]
    astCount = astObj.sum()
    astSignal = astCount-(avgSky*((astLength+4)**2))
    astMag = -2.5*math.log10(astSignal)
    print("The asteroid's instrumental mag is ", astMag)

def calcSky(x, y, inner, outer):
    return
    
    
fractionalPhotometry("sampleimage.fits")
