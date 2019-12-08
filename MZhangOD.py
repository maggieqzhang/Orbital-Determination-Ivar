# Method of Gauss

from Functions import *
import numpy as np
from math import *

# change colons to spaces in file
with open("MZhangInput.txt", "r+") as data:
    txt = data.read()
    txt = txt.replace(':', " ")

file = open("MZhangInput.txt", "w")
file.write(txt)
file.close()

# load data from input file
data = np.array(np.loadtxt("MZhangInput.txt"))
        
# initial values from input file
dates = data[:,:6]
RA = data[:,6:9]
dec = data[:,9:12]
Rx = data[:,12:13]
Ry = data[:,13:14]
Rz = data[:,14:]
k = 0.01720209895

# conversion of input data to proper format
dates = gaussian(dates) # in gaussian days
RA = raToRad(RA) # ra in radians
dec = decToRad(dec) # dec in radians

# Turning R vectors into a compiled matrix
R = np.array([[]])
for i in range(len(RA)):
    temp = np.array([[Rx[i], Ry[i], Rz[i]]])
    R = np.append(R, temp)
R = np.reshape(R, ((len(RA), 3)))

# uncertainty on RA and Dec are manually inputted 
uncertaintyRa = np.array([0.7666432826405849,2.965564005970234,2.0272810546312785,1.3164310034687199])
uncertaintyRa = (uncertaintyRa/3600)
uncertaintyDec = np.array([0.6472287404975303,2.9045272981816663,1.7159101550688076, 1.006005772642148])
uncertaintyDec = (uncertaintyDec/3600)

# specifying which observations to use
for i in range(1, data.shape[0]-1):
    obs1 = 0 
    obs2 = i
    obs3 = data.shape[0]-1

    # finding sun vectors
    R1 = np.array([Rx[obs1][0], Ry[obs1][0], Rz[obs1][0]])
    R2 = np.array([Rx[obs2][0], Ry[obs2][0], Rz[obs2][0]])
    R3 = np.array([Rx[obs3][0], Ry[obs3][0], Rz[obs3][0]])

    # finding rhohat vectors
    rhohat1 = rhohat(RA[obs1], dec[obs1])
    rhohat2 = rhohat(RA[obs2], dec[obs2])
    rhohat3 = rhohat(RA[obs3], dec[obs3])

    # finding D constants
    D0 = np.dot(rhohat1, np.cross(rhohat2, rhohat3))
    D11 = np.dot(np.cross(R1, rhohat2), rhohat3)
    D12 = np.dot(np.cross(R2, rhohat2), rhohat3)
    D13 = np.dot(np.cross(R3, rhohat2), rhohat3)
    D21 = np.dot(np.cross(rhohat1, R1), rhohat3)
    D22 = np.dot(np.cross(rhohat1, R2), rhohat3)
    D23 = np.dot(np.cross(rhohat1, R3), rhohat3)
    D31 = np.dot(rhohat1, np.cross(rhohat2, R1))
    D32 = np.dot(rhohat1, np.cross(rhohat2, R2))
    D33 = np.dot(rhohat1, np.cross(rhohat2, R3))
    
    # find possible r2 values
    r2 = scalarEqOfLagrange(dates[obs1], dates[obs2], dates[obs3], rhohat2, R2, D21, D22, D23, D0)
    print("Possible Starting r2 values :", r2)
    rootChosen = int(input("Which root would you like to use?"))
    print("")
    root = r2[rootChosen]
    mu = 1
    tau1 = dates[obs1] - dates[obs2]
    tau3 = dates[obs3] - dates[obs2]
    u2 = mu/root**3
    f1 = 1 - (mu*tau1**2 / (2*root**3))
    f3 = 1 - (mu*tau3**2 / (2*root**3))
    g1 = tau1 - (mu*tau1**3/(6*root**3))
    g3 = tau3 - (mu*tau3**3/(6*root**3))
    c1 = g3/(f1*g3-f3*g1)
    c2 = -1
    c3 = -1*g1/(f1*g3-f3*g1)
    rho1 = (c1*D11+c2*D12+c3*D13)/(c1*D0)
    rho2 = (c1*D21+c2*D22+c3*D23)/(c2*D0)
    rho3 = (c1*D31+c2*D32+c3*D33)/(c3*D0)
    r1 = rho1 * rhohat1 - R1
    r2 = rho2 * rhohat2 - R2
    r3 = rho3 * rhohat3 - R3
    d1 = -1*f3/(f1*g3-f3*g1)
    d3 = f1/(f1*g3-f3*g1)
    r2dot = d1*r1 + d3*r3

    # prints original 
    print("Original r2 and r2dot vectors")
    print("r2 vector :", r2.real, "AU")
    print("r2 dot vector :", k*r2dot.real, "AU/day")
    print("")

    # iterates with Taylor Series
    finished = iterateTaylor(root.real, r2.real, r2dot.real, dates[obs1], dates[obs2], dates[obs3], rhohat1, rhohat2, rhohat3, R1, R2, R3, D0, D11, D12, D13, D21, D22, D23, D31, D32, D33)
    print("After Method of Gauss r2 and r2dot")
    print("r2 vector :", eqToEc(finished[0]), "AU")
    print("r2dot vector :", k*eqToEc(finished[1]), "AU/day")
    print("")
    print("Range to asteroid :", findMag(finished[0]), "AU")
    print("")
    oe = orbitalElements(finished[0], finished[1], dates[obs2]/k)
    # Orbital
    print("Orbital Elements with MoG r2 and r2dot")
    print("a =",oe[1], "AU")
    print("e =",oe[0])
    print("i = " + str(degrees(oe[2])) + "°")
    print("Ω = "+ str(degrees(oe[3]))+"°")
    print("ω = " + str(degrees(oe[4])%360) + "°")
    print("M =", str(degrees(oe[5])) + "°")
    pairOrig = generateEphem(oe[0], oe[1], oe[2], oe[3], oe[4], oe[5], oe[6], oe[7], R2)

    # differential correction
    diffDates = np.copy(dates)
    corrected = differentialCorrection(RA, dec, diffDates, oe[0], oe[1], oe[2], oe[3], oe[4], oe[5], oe[6], R, finished[0], finished[1])
    orbital = orbitalElements(corrected[:3] + finished[0], corrected[3:]+finished[1], dates[obs2]/k)
    print("")
    print("Orbital Elements after Differential Correction")
    print("a =",orbital[1], "AU")
    print("e =",orbital[0])
    print("i = " + str(degrees(orbital[2])) + "°")
    print("Ω = "+ str(degrees(orbital[3]))+"°")
    print("ω = " + str(degrees(orbital[4])%360) + "°")
    print("M =", str(degrees(orbital[5])) + "°")
    print("")

# Monte Carlo
run = str(input("Run Monte Carlo 100,000 times? 'y' or 'n'"))
if run == "y":
    monteDates = np.copy(dates)
    monteCarlo(RA, dec, uncertaintyRa, uncertaintyDec, monteDates, R, 1)
    
