# OD LSPR code
# Maggie Zhang

import numpy as np
import astropy
from astropy.io import fits
from math import *
import matplotlib.pyplot as plt

sample = np.array(np.loadtxt("july7set1.txt"))

astX = float(input("X coord of the asteroid"))
astY = float(input("Y coord of the asteroid"))
flatten = str(input("Do you want to flatten? 'y' or 'n'"))
    
N = 12
x = sample[:,:1]
y = sample[:,1:2]
ra = sample[:,2:3]
dec = sample[:,3:]

# round to a number of sig figs
def roundSig(x, sig):
    return round(x, sig-int(floor(log10(x)))-1)

# xy sums
sumX = 0
sumY = 0
sumXSquared = 0
sumxy = 0
sumYSquared = 0

for i in range(12):
    sumX = sumX + x[i]
    sumY = sumY + y[i]
    sumXSquared = sumXSquared + (x[i])**2
    sumxy = sumxy + x[i]*y[i]
    sumYSquared = sumYSquared + (y[i])**2

# ra dec sums
sumDec = 0
sumDecX = 0
sumDecY = 0
sumRA = 0
sumRAX = 0
sumRAY = 0

for i in range(12):
    sumDec = sumDec + dec[i]
    sumRA = sumRA + ra[i]
    sumDecX = sumDecX + dec[i]*x[i]
    sumDecY = sumDecY + dec[i]*y[i]
    sumRAX = sumRAX + ra[i]*x[i]
    sumRAY = sumRAY + ra[i]*y[i]

# RA matrices
raSums = np.array([sumRA, sumRAX, sumRAY]).reshape((3,1))
A = np.array([[N, sumX, sumY],[sumX, sumXSquared, sumxy],[sumY, sumxy, sumYSquared]])
A = np.linalg.inv(A)
raConstants = np.dot(A, raSums)
b1 = raConstants[0][0]
a11 = raConstants[1][0]
a12 = raConstants[2][0]

# Dec matrices
decSums = np.array([sumDec, sumDecX, sumDecY]).reshape((3,1))
decConstants = np.dot(A, decSums)
b2 = decConstants[0][0]
a21 = decConstants[1][0]
a22 = decConstants[2][0]


#uncertainty
fitRa = []
fitDec = []
for i in range(12):
    newRa = b1 + a11*x[i] + a12*y[i]
    newDec = b2 + a21*x[i] + a22*y[i]
    fitRa.append(newRa)
    fitDec.append(newDec)

sigmaRa = 0
sigmaDec = 0
squaredRa = 0
squaredDec = 0
for i in range(12):
    squaredRa =  squaredRa + (ra[i]-fitRa[i])**2
    squaredDec = squaredDec + (dec[i]-fitDec[i])**2
sigmaRa = sqrt(squaredRa/(N-3))
sigmaDec = sqrt(squaredDec/(N-3))


def formatAngle(angle):
    deg = int(angle)
    minute = (angle - deg)*60
    degmin = int(minute)
    second = (minute-degmin)*60
    return str(deg) + ":" + str(degmin) + ":" + str(second)
    
# asteroid
astRA = b1 + a11*astX + a12*astY
astDec = b2 + a21*astX + a22*astY

if flatten == "n":
#printing plate constants
    print("Plate Contants")
    print("b1: " + str(b1) + " deg")
    print("b2: " + str(b2) + " deg")
    print("a11: " + str(a11) + " deg/pix")
    print("a12: " + str(a12) + " deg/pix")
    print("a21: " + str(a21) + " deg/pix")
    print("a22: " + str(a22) + " deg/pix")
    print()
    print("Uncertainty")
    print("RA: " + str(sigmaRa * 3600) + " arcsec")
    print("Dec: " + str(sigmaDec * 3600) + " arcsec")
    print()
    print("Astronomy for (x,y) = (" + str(astX) + ", " + str(astY) + ")")
    print("RA= " + formatAngle(astRA*24/360))
    print("Dec= +" + formatAngle(astDec))
    
##hdulist = fits.open('REAL IVAR.00000216.ENTERED COORDINATES.REDUCED.FIT', mode = 'update')
##hdu = hdulist[0].header
##hdu.set('CTYPE1', 'RA---TAN')
##hdu.set('CRPIX1', float(sumX/N))
##hdu.set('CRVAL1', b1)
##hdu.set('CTYPE2', 'DEC--TAN')
##hdu.set('CRPIX2', float(sumY/N))
##hdu.set('CRVAL2', b2)
##hdu.set('CD1_1', a11)
##hdu.set('CD1_2', a12)
##hdu.set('CD2_1', a21)
##hdu.set('CD2_2', a22)
##hdu.set('RADECSYS', 'FK5')
##hdu.set('EQUINOX', 2000.0)
##hdulist.close()

# Flattening
AConstant = radians(sumRA/12)
D = radians(sumDec/12)
L = 3911/0.024

def flattenDec(ra, dec, y):
    ra = radians(ra)
    dec = radians(dec)
    H=sin(dec)*sin(D)+cos(dec)*cos(D)*cos(ra-AConstant)
    flatDec = ((sin(dec)*cos(D)-cos(dec)*sin(D)*cos(ra-AConstant))/H)-(y/L)
    return flatDec

def flattenRA(ra, dec, x):
    ra = radians(ra)
    dec = radians(dec)
    H=sin(dec)*sin(D)+cos(dec)*cos(D)*cos(ra-AConstant)
    flatRa = (cos(dec)*sin(ra-AConstant)/H)-(x/L)
    return flatRa

flatDecs = []
flatRas = []
for i in range(12):
    flatRas.append(flattenRA(ra[i], dec[i], x[i]))
    flatDecs.append(flattenDec(ra[i], dec[i], y[i]))

sumFDec = 0
sumFDecX = 0
sumFDecY = 0
sumFRA = 0
sumFRAX = 0
sumFRAY = 0

for i in range(12):
    sumFDec = sumFDec + flatDecs[i]
    sumFRA = sumFRA + flatRas[i]
    sumFDecX = sumFDecX + flatDecs[i]*x[i]
    sumFDecY = sumFDecY + flatDecs[i]*y[i]
    sumFRAX = sumFRAX + flatRas[i]*x[i]
    sumFRAY = sumFRAY + flatRas[i]*y[i]

# RA matrices
raFlatSums = np.array([sumFRA, sumFRAX, sumFRAY]).reshape((3,1))
raConstants = np.dot(A, raFlatSums)
fb1 = raConstants[0][0]
fa11 = raConstants[1][0]
fa12 = raConstants[2][0]

decFlatSums = np.array([sumFDec, sumFDecX, sumFDecY]).reshape((3,1))
decConstants = np.dot(A, decFlatSums)
fb2 = decConstants[0][0]
fa21 = decConstants[1][0]
fa22 = decConstants[2][0]


flatastRA = fb1 + fa11*astX + fa12*astY + astX/L
flatastDec = fb2 + fa21*astX + fa22*astY + astY/L

def fromFlatRa(ra, dec):
    delta = cos(D)-dec*sin(D)
    RA = AConstant + atan(ra/delta)
    return RA

def fromFlatDec(ra, dec):
    delta = cos(D)-dec*sin(D)
    r = sqrt(ra**2 + delta**2)
    Dec = atan((sin(D)+dec*cos(D))/r)
    return Dec

fitflatRa = []
fitflatDec = []
for i in range(12):
    newfRa = fb1 + fa11*x[i] + fa12*y[i] + x[i]/L
    newfDec = fb2 + fa21*x[i] + fa22*y[i] + y[i]/L
    newRa = fromFlatRa(newfRa, newfDec)
    newDec = fromFlatDec(newfRa, newfDec)
    fitflatRa.append(degrees(newRa))
    fitflatDec.append(degrees(newDec))

fAstRa = fromFlatRa(flatastRA, flatastDec)
fAstDec = fromFlatDec(flatastRA, flatastDec)


fsigmaRa = 0
fsigmaDec = 0
fsquaredRa = 0
fsquaredDec = 0
for i in range(12):
    fsquaredRa =  fsquaredRa + (ra[i]-((fitflatRa[i])))**2
    fsquaredDec = fsquaredDec + (dec[i]-((fitflatDec[i])))**2
fsigmaRa = sqrt(fsquaredRa/(N-3))
fsigmaDec = sqrt(fsquaredDec/(N-3))

if flatten == "y":
    print()
    print("Flattened Plate Contants")
    print("b1: " + str(fb1))
    print("b2: " + str(fb2))
    print("a11: " + str(fa11))
    print("a12: " + str(fa12))
    print("a21: " + str(fa21))
    print("a22: " + str(fa22))
    print()
    print("Flattened Astrometry for (x,y) = (" + str(astX) + ", " + str(astY) + ")")
    print("RA = " + str(formatAngle(degrees(fAstRa/15))))
    print("Dec = +" + str(formatAngle(degrees(fAstDec))))
    print()
    print("Flattened Uncertainty")
    print("RA: " + str(roundSig(fsigmaRa * 3600, 1)) + " arcsec")
    print("Dec: " + str(roundSig(fsigmaDec * 3600, 1)) + " arcsec")

