# Math/Physics PSet 5
# Maggie Zhang

from math import *
import numpy as np

# functions
def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return asin(sine)

    if cosine < 0 and sine > 0: #2
        return acos(cosine)

    if cosine < 0 and sine < 0: #3
        return pi - asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*pi + asin(sine)

def eqToEc(a):
    epsilon = radians(-23.4352)
    ax = a[0]
    ay = a[1]
    az = a[2]
    newx = ax
    newy = ay*cos(epsilon) - az*sin(epsilon)
    newz = ay*sin(epsilon) + az*cos(epsilon)
    return np.array([newx, newy, newz]) 

def solvekep(M, e):
    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    Mtrue = M
    while abs(Mguess - Mtrue) > 1e-004:
        Mguess = Eguess - e*sin(Eguess)
        Eguess = Eguess - (Mguess-Mtrue) / (1-e*cos(Eguess))
    return Eguess

def formatAngle(angle):
    deg = int(angle)
    minute = (angle - deg)*60
    degmin = int(minute)
    second = (minute-degmin)*60
    return str(deg) + ":" + str(degmin) + ":" + str(second)

def generateEphem(e, a, i, omega, w, M0, t0, t):
    epsilon = radians(23.4352)
    rootmu = 0.01720209895
    # mean motion
    n = rootmu*sqrt(1/a**3)

    # mean anomaly
    M = M0 + n*(t-t0)

    # eccentric anomaly 
    E = solvekep(M, e)

    # x and y coordinates
    x = a*(cos(E)-e)
    y = a * sin(E) * sqrt(1-e**2)
    r = sqrt(x**2 + y**2)

    # true anomaly
    f = findQuadrant(y/r, x/r)

    # asteroid ecliptic coordinates
    rX = r*(cos(f+w)*cos(omega) - cos(i)*sin(f+w)*sin(omega))
    rY = r*(cos(i)*cos(omega)*sin(f+w) + cos(f+w)*sin(omega))
    rZ = r*(sin(i) * sin(f+w))
    print(rX, rY, rZ)

    # sun ecliptic coords for R from jpl horizons
    R = np.array([-1.761362304754478E-01, 9.186910039895203E-01 ,3.982195186381176E-01])
    R = eqToEc(R)
    RX = R[0]
    RY = R[1]
    RZ = R[2]

    # components of earth to asteroid vector rho
    rhoX = rX + RX
    rhoY = rY + RY
    rhoZ = rZ + RZ

    # components of rho equatorial
    rhoXEq = rhoX
    rhoYEq = rhoY*cos(epsilon) - rhoZ*sin(epsilon)
    rhoZEq = rhoY*sin(epsilon) + rhoZ*cos(epsilon)
    rho = sqrt(rhoX**2 + rhoY**2 + rhoZ**2)

    # find ra and dec
    dec = asin(rhoZEq/rho)
    sinRA = rhoYEq / (rho * cos(dec))
    cosRA = rhoXEq /(rho*cos(dec))
    ra = findQuadrant(sinRA, cosRA)

    print("Dec =", formatAngle(degrees(dec)))
    print("RA =", formatAngle(degrees(ra)/15))

e =  0.39378153749170036
a = 1.8534126044942794 #AU
i = radians( 8.406512659898615)
omega = radians(133.17295119651467)
w = radians(167.76478487779573)
M = radians(353.5214322928663)
rootmu = 0.01720209895
t0 = 2458321.75#JHD
t = 2458301.648330 #JHD
generateEphem(e, a, i, omega, w, M, t0, t)



























