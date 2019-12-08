# Monte Carlo Simulation

import numpy as np
from math import *
import Functions
import matplotlib.pyplot as plt
from matplotlib import colors

# lists of a, e, i, omega, w, M values from running the Method of Gauss and baby OD
a = np.array([])
ec = np.array([])
inclination = np.array([])
omega = np.array([])
w = np.array([])
M = np.array([])
count = 0
        
def runMonte(inputFile, obsNum):
    global a, ec, inclination, omega, w, M, count #, a2, ec2, i2, omega2, w2, M2

    count+=1
    data = np.array(np.loadtxt(inputFile))
    dates = data[:,:1]
    RA = data[:,1:2]
    dec = data[:,2:3]
    Rx = data[:,3:4]
    Ry = data[:,4:5]
    Rz = data[:,5:]
    k = 0.01720209895

    R = np.array([[]])
    for i in range(len(RA)):
        temp = np.array([[Rx[i], Ry[i], Rz[i]]])
        R = np.append(R, temp)
    R = np.reshape(R, ((len(RA), 3)))

    # specifying which observations to use
    obs1 = 0 
    obs2 = 1
    obs3 = 2

    # finding sun vectors
    R1 = np.array([Rx[obs1][0], Ry[obs1][0], Rz[obs1][0]])
    R2 = np.array([Rx[obs2][0], Ry[obs2][0], Rz[obs2][0]])
    R3 = np.array([Rx[obs3][0], Ry[obs3][0], Rz[obs3][0]])

    # finding rhohat vectors
    rhohat1 = Functions.rhohat(RA[obs1], dec[obs1])
    rhohat2 = Functions.rhohat(RA[obs2], dec[obs2])
    rhohat3 = Functions.rhohat(RA[obs3], dec[obs3])

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
    r2 = Functions.scalarEqOfLagrange(dates[obs1], dates[obs2], dates[obs3], rhohat2, R2, D21, D22, D23, D0)
    # run through each root from the Scalar Equation of Lagrange 
    for root in r2:
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
        finished = Functions.iterateTaylor(root.real, r2.real, r2dot.real, dates[obs1], dates[obs2], dates[obs3], rhohat1, rhohat2, rhohat3, R1, R2, R3, D0, D11, D12, D13, D21, D22, D23, D31, D32, D33)
        if finished == None:
            return
        oe = Functions.orbitalElements(finished[0], finished[1], dates[obs2]/k)
        a = np.append(a, oe[1])
        ec = np.append(ec, oe[0])
        inclination = np.append(inclination, degrees(oe[2]))
        omega = np.append(omega, degrees(oe[3]))
        w = np.append(w, (degrees(oe[4]))%360)
        M = np.append(M, degrees(oe[5]))
        # check to see when count is the same as the number of iterations of the method is run
        if count == 10000: # 10,000 runs through monte carlo
            count = 1
            plotMonte(a, "Semi-Major Axis ", "(AU)")
            plotMonte(ec, "Eccentricity", "")
            plotMonte(inclination, "Inclination ", "(°)")
            plotMonte(omega, "Longitude of the Ascending Node ", "(°)")
            plotMonte(w, "Argument of Perihelion ", "(°)")
            plotMonte(M, "Mean Anomaly ", "(°)")

# plots histograms of each orbital element
def plotMonte(vals, title, unit):
    print(title, "mean", np.sum(vals)/len(vals))
    print(title, "standard deviation", np.std(vals))
    Na, binsa, patchesa = plt.hist(vals, bins = 50)
    fraca = Na/Na.max()
    norm = colors.Normalize(fraca.min(), fraca.max())
    for fraction, patch in zip(fraca, patchesa):
        color = plt.cm.viridis(norm(fraction))
        patch.set_facecolor(color)
    plt.xlabel(title + unit, fontsize=10, fontname="Times New Roman")
    plt.ylabel("Frequency (N)", fontsize=10, fontname="Times New Roman")
    plt.title("Monte Carlo Distribution of " + title, fontname="Times New Roman")
    plt.show()


