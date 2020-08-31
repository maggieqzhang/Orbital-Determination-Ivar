# OD Centroid

#Centroiding with FITS Image
#Maggie Zhang
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

def centroid(image, data):
    im = fits.getdata(image)
    x = data[:,:1]
    y = data[:,1:2]
    lengthX = data[:,2:3]
    lengthY = data[:,3:4]
    ra = data[:,4:5]
    ra = 15*ra
    dec = data[:,5:]
    for i in range(12):
        halfX = int(lengthX[i]/2)
        halfY = int(lengthY[i]/2)
        star = 0
        Y = int(y[i])
        X = int(x[i])
        star = im[Y-halfY:Y+halfY+1, X-halfX:X+halfX+1]
        if lengthX[i]%2 == 0:
            if lengthY[i]%2==0:
                star = im[Y-halfY:Y+halfY, X-halfX:X+halfX]
            else:
                star = im[Y-halfY:Y+halfY+1, X-halfX:X+halfX]
        else:
            if lengthY[i]%2==0:
                star = im[Y-halfY:Y+halfY, X-halfX:X+halfX+1]
            else:
                star = im[Y-halfY:Y+halfY+1, X-halfX:X+halfX+1]
        totalPixelCount = star.sum()
        weightedX = 0
        weightedY = 0
        for index in range(star.shape[1]):
            weightedX = weightedX + star[:,index:(index+1)].sum()*index
        for index in range(star.shape[0]):
            weightedY = weightedY + star[index].sum()*index
        centroidX = weightedX/totalPixelCount+X-1
        centroidY = weightedY/totalPixelCount+Y-1
        centroidNumbers.write(str(centroidX))
        centroidNumbers.write(" ")
        centroidNumbers.write(str(centroidY))
        centroidNumbers.write(" ")
        centroidNumbers.write(str(ra[i][0]))
        centroidNumbers.write(" ")
        centroidNumbers.write(str(dec[i][0]))
        centroidNumbers.write("\n")
        
##data = np.array(
##[[371,489,12,11,14.79779857,5.72922583],
##[239,608,13,11,14.80080915,5.68772583],
##[527,684,11,11,14.79415663,5.6629625],
##[221,788,7,7,14.80116009,5.62559111],
##[974,606,10,9,14.78385676,5.69131],
##[921,479,11,7,14.787742,5.73289972],
##[612,437,9,9,14.79224335,5.74817528],
##[657,486,5,5,14.79120824,5.73119306],
##[380,623,6,5,14.79753728,5.68283833],
##[442,419,7,5,14.79619826,5.75337139],
##[488,503,4,4,14.79509294,5.72484722],
##[499,230,5,5,14.79491948,5.81868528]])

centroids = np.array(np.loadtxt("centroids.txt."))
centroidNumbers = open("July18.txt", "w")
centroid("IVAR LAST DITCH.00000083.ENTERED COORDINATES.REDUCED.FIT", centroids)
centroidNumbers.close()

