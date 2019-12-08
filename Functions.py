# Method of Gauss and Monte Carlo functions

from math import *
import numpy as np
from montecarlocode import *
import time

def gaussian(ut):
    gauss = []
    k = 0.01720209895
    for i in range(ut.shape[0]):
        year = ut[i][0]
        month = ut[i][1]
        day = ut[i][2]
        hour = ut[i][3]
        minute = ut[i][4]
        second = ut[i][5]
        # j0 formula
        j0 = 367*year - int(7*(year + int((month+9)/12))/4) + int(275*month/9) + day + 1721013.5
        decimalHour = hour + minute/60 + second/3600
        jhd = j0 + decimalHour/24
        gDays = jhd*k
        gauss.append(gDays)
    return gauss

def raToRad(ra):
    RA = []
    for i in range(ra.shape[0]): # shape - rows by cols
        decimal = ra[i][0]
        decimal += ra[i][1]/60
        decimal += ra[i][2]/3600
        rad = radians(decimal * 15)
        RA.append(rad)
    return RA

def decToRad(dec):
    DEC = []
    for i in range(dec.shape[0]):
        decimal = dec[i][0]
        decimal += dec[i][1]/60
        decimal += dec[i][2]/3600
        rad = radians(decimal)
        DEC.append(rad)
    return DEC

def rhohat(ra, dec):
    x = cos(ra) * cos(dec)
    y = sin(ra) * cos(dec)
    z =  sin(dec)
    rhoHat = np.array([x, y, z])
    return rhoHat

def findMag(vec):
    mag = vec[0]**2 + vec[1]**2 + vec[2]**2
    mag = sqrt(mag)
    return mag

def scalarEqOfLagrange(tau1, tau2, tau3, rhohat2, R2, D21, D22, D23, D0):
    k = 0.01720209895
    mu = 1

    Tau1 = (tau1-tau2)
    Tau3 = (tau3-tau2)
    tau = tau3-tau1
    
    A1 = Tau3/tau
    B1 = A1*(tau**2-Tau3**2)/6
    A3 = -1*Tau1/tau
    B3 = A3*(tau**2-Tau1**2)/6
    
    A = (A1*D21 - D22 + A3*D23) / (-1*D0)
    B = (B1*D21 + B3*D23)/(-1*D0)
    E = -1*2*(np.dot(rhohat2, R2))
    F = (findMag(R2))**2
    
    a = -1*(A**2 + A * E + F)
    b = -1*mu * (2*A*B + B*E)
    c = -1*mu**2 * B**2
    
    coeff = [1, 0, a, 0, 0, b, 0, 0, c]
    roots = np.roots(coeff)
    roots = roots[np.isreal(roots)]
    finalroots = []
    
    # checking to see if rho value and r are negative
    for r in roots:
        rhor = A + mu*B/r**3
        if r > 0 and rhor > 0:
            finalroots.append(r)
    return finalroots

# finding Eccentric Anomaly for closed form iterations
def solveKep(nt, e, E2):
    Eguess0 = nt
    Xold = 0
    X = Eguess0
    while abs(Xold-X) > 1e-004:
        Xold = X
        X = X - (X - e*cos(E2)*sin(X) + e*sin(E2)*(1-cos(X))-Eguess0)/(1-e*cos(E2)*cos(X)+e*sin(E2)*sin(X))
    return X

# finding Eccentric Anomaly for ephemeris generation using Mean Anomaly
def solvekep1(M, e):
    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    Mtrue = M
    while abs(Mguess - Mtrue) > 1e-004:
        Mguess = Eguess - e*sin(Eguess)
        Eguess = Eguess - (Mguess-Mtrue) / (1-e*cos(Eguess))
    return Eguess

def eqToEc(a):
    epsilon = radians(-23.4352)
    ax = a[0]
    ay = a[1]
    az = a[2]
    newx = ax
    newy = ay*cos(epsilon) - az*sin(epsilon)
    newz = ay*sin(epsilon) + az*cos(epsilon)
    return np.array([newx, newy, newz]) 

# closed form interation in the Method of Gauss
def iterateClosed(r2Original, r20, r2dot0, tau1, tau2, tau3, rhohat1, rhohat2, rhohat3, R1, R2, R3, D0, D11, D12, D13, D21, D22, D23, D31, D32, D33):
    # constants
    k = 0.01720209895
    c = 173.145 # speed of light in au/day
    mu = 1
    count = 0

    # time calculations
    t1 = tau1/k # Original t in JHD
    t2 = tau2/k # Original t in JHD
    t3 = tau3/k # Original t in JHD
    T1 = t1 # changing t in JHD
    T2 = t2 # changing t in JHD
    T3 = t3 # changing t in JHD
    Tau1 = k*(T1-T2)
    Tau3 = k*(T3-T2)
    Tau = k*(T3-T1)

    # finding magnitudes and defining the final r2 and r2dot to return
    r2 = r20
    rmagprevious = r2Original
    rmagnew = findMag(r20)
    r2dot = r2dot0
    finalr2 = r2
    finalr2dot = r2dot

    #begin iteration loop
    while abs(rmagprevious - rmagnew) > 10**(-12):
        count+= 1
        rmagprevious = rmagnew
        # constant orbital elemnts to find f and g
        a = ((2/findMag(eqToEc(r2))) - np.dot(eqToEc(r2dot), eqToEc(r2dot)))**(-1)
        e = sqrt((1-np.dot(np.cross(eqToEc(r2),eqToEc(r2dot)), np.cross(eqToEc(r2),eqToEc(r2dot)))/a))
        n = k * sqrt(mu/a**3) # use k = 0.017 because time is in JHD
        M1 = n * (T1-T2)
        M3 = n * (T3-T2)
        E2 = acos((1-findMag(eqToEc(r2))/a)/e)
        E1 = solveKep(M1, e, E2)
        E3 = solveKep(M3, e, E2)
        
        nGauss = sqrt(mu/a**3) # use mu = 1 because everything is in Gaussian Days
        f1 = 1 - a*(1-cos(E1))/rmagprevious
        f3 = 1 - a*(1-cos(E3))/rmagprevious
        g1 = Tau1 + (sin(E1)-(E1))/nGauss
        g3 = Tau3 + (sin(E3)-(E3))/nGauss
        
        c1 = g3/(f1*g3-f3*g1)
        c2 = -1
        c3 = -1*g1/(f1*g3-f3*g1)
        
        rho1 = (c1*D11+c2*D12+c3*D13)/(c1*D0)
        rho2 = (c1*D21+c2*D22+c3*D23)/(c2*D0)
        rho3 = (c1*D31+c2*D32+c3*D33)/(c3*D0)
        
        r1 = rho1 * rhohat1 - R1
        r2 = rho2 * rhohat2 - R2
        r3 = rho3 * rhohat3 - R3
        
        rmagnew = findMag(r2)
        finalr2 = r2
        
        d1 = -1*f3/(f1*g3-f3*g1)
        d3 = f1/(f1*g3-f3*g1)
        r2dot = d1*r1 + d3*r3
        finalr2dot = r2dot
        
        # time light correction 
        T1 = t1 - rho1/c
        T2 = t2 - rho2/c
        T3 = t3 - rho3/c
        Tau1 = k*(T1-T2)
        Tau3 = k*(T3-T2)
    # if the series does not converge after 1000 iterations, exit
    if count == 1000:
        return None
    return finalr2, finalr2dot

# Taylor series interation in the Method of Gauss
def iterateTaylor(r2Orig, r20, r2dot0, tau1, tau2, tau3, rhohat1, rhohat2, rhohat3, R1, R2, R3, D0, D11, D12, D13, D21, D22, D23, D31, D32, D33):
    # constants
    k = 0.01720209895
    c = 173.145 # speed of light in au/day
    mu = 1
    count = 0

    # time calculations
    t1 = tau1/k # Original t in JHD
    t2 = tau2/k # Original t in JHD
    t3 = tau3/k # Original t in JHD
    T1 = t1 # changing t in JHD
    T2 = t2 # changing t in JHD
    T3 = t3 # changing t in JHD
    Tau1 = k*(T1-T2)
    Tau3 = k*(T3-T2)
    Tau = k*(T3-T1)

    # finding magnitudes and defining the final r2 and r2dot to return
    r2 = r20
    rmagprevious = r2Orig
    rmagnew = findMag(r20)
    r2dot = r2dot0
    finalr2 = r2
    finalr2dot = r2dot

    #begin iteration loop
    while abs(rmagprevious - rmagnew) > 10**(-12) and count < 1000:
        count +=1
        rmagprevious = rmagnew
        
        u = mu/rmagprevious**3 # use mu = 1 because times in f and g series are in Gaussian Days
        q = (np.dot(r2dot, r2dot))/rmagprevious**2 - u
        z = (np.dot(r2, r2dot))/rmagprevious**2
        f1 = 1- (u*Tau1**2)/2 + (u*z*Tau1**3)/2 + ((3*u*q-15*u*z**2+u**2)/24)*Tau1**4
        f3 = 1- (u*Tau3**2)/2 + (u*z*Tau3**3)/2 + ((3*u*q-15*u*z**2+u**2)/24)*Tau3**4
        g1 = Tau1 - (u*Tau1**3)/6 + (u*z*Tau1**4)/4
        g3 = Tau3 - (u*Tau3**3)/6 + (u*z*Tau3**4)/4
        
        c1 = g3/(f1*g3-f3*g1)
        c2 = -1
        c3 = -1*g1/(f1*g3-f3*g1)
        
        rho1 = (c1*D11+c2*D12+c3*D13)/(c1*D0)
        rho2 = (c1*D21+c2*D22+c3*D23)/(c2*D0)
        rho3 = (c1*D31+c2*D32+c3*D33)/(c3*D0)
        
        r1 = rho1 * rhohat1 - R1
        r2 = rho2 * rhohat2 - R2
        r3 = rho3 * rhohat3 - R3
        
        rmagnew = findMag(r2)
        finalr2 = r2
        
        d1 = -1*f3/(f1*g3-f3*g1)
        d3 = f1/(f1*g3-f3*g1)
        r2dot = d1*r1 + d3*r3
        finalr2dot = r2dot

        # time light correction 
        T1 = t1 - rho1/c
        T2 = t2 - rho2/c
        T3 = t3 - rho3/c
        Tau1 = k*(T1-T2)
        Tau3 = k*(T3-T2)
    # if the series does not converge after 1000 iterations, exit
    if count == 1000:
        return None
    return finalr2, finalr2dot    

def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return asin(sine)

    if cosine < 0 and sine > 0: #2
        return acos(cosine)

    if cosine < 0 and sine < 0: #3
        return pi - asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*pi + asin(sine)

def formatAngle(angle):
    deg = int(angle)
    minute = (angle - deg)*60
    degmin = int(minute)
    second = (minute-degmin)*60
    return str(deg) + ":" + str(degmin) + ":" + str(second)
    
def orbitalElements(r2, r2dot, t2):
    rootmu = 0.01720209895
    r2 = eqToEc(r2)
    r2dot = eqToEc(r2dot)

    r2x = r2[0]
    r2y = r2[1]
    r2z = r2[2]

    r2dotx = r2dot[0]
    r2doty = r2dot[1]
    r2dotz = r2dot[2]

    r = sqrt(r2x**2 + r2y**2 + r2z**2)
    rdot = sqrt(r2dotx**2 + r2doty**2 + r2dotz**2)

    # a - semi-major axis
    a = (2/r - np.dot(r2dot, r2dot))**-1

    # e - eccentricity
    e = sqrt((1-np.dot(np.cross(r2, r2dot), np.cross(r2, r2dot))/a))

    # i - inclination
    h = rootmu * np.cross(r2, r2dot)
    hx = h[0] 
    hy = h[1] 
    hz = h[2]
    hmag = sqrt(hx**2 + hy**2 + hz**2)
    i = acos(hz/hmag)

    # omega
    sinOmega = hx/(hmag*sin(i))
    cosOmega = -hy/(hmag*sin(i))
    omega = findQuadrant(sinOmega, cosOmega)

    # w
    sinwf = r2z/(r*sin(i))
    coswf = 1/cos(omega)*((r2x/r) + cos(i)*sinwf*sin(omega))
    wf = findQuadrant(sinwf, coswf)
    sinf = (np.dot(r2,r2dot))*sqrt(a*(1-e**2))/(e*r)
    cosf = (a*(1-e**2)/r - 1)/e
    f = findQuadrant(sinf,cosf)
    w = wf-f

    # M0 - mean anomaly
    E2 = acos((1-r/a)/e)
    M2 = E2 - e*sin(E2)
    if (M2 > pi and f < pi) or (M2 < pi and f > pi):
        M2 = 2* pi - M2
    k = 0.01720209895
    mu = 1
    n = k * sqrt(mu/a**3)
    t = gaussian(np.array([[2018, 7, 22, 6, 0, 0]]))[0]
    M = M2 + n*(t/k- t2)
    return e, a, i, omega, w, M, (gaussian(np.array([[2018, 7, 22, 6, 0, 0]]))[0])/k, t2

def generateEphem(e, a, i, omega, w, M0, t0, t, R):
    epsilon = radians(23.4352)
    rootmu = 0.01720209895
    R = eqToEc(R)
    # mean motion
    n = rootmu*sqrt(1/a**3)

    # mean anomaly
    M = M0 + n*(t-t0)

    # eccentric anomaly 
    E = solvekep1(M, e)

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

    # sun ecliptic coords for R from jpl horizons
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
    return ra, dec

def differentialCorrection(ra, dec, t, e, a, i, omega, w, M0, t0, R, r2, r2dot):
    delta = 10**-3
    for index in range(len(t)):
        time = t[index] /0.01720209895
        t[index] = time
    calculated = np.array([])
    radec = np.array([])
    for index in range(len(ra)):
        pair = generateEphem(e, a, i, omega, w, M0, t0, t[index], R[index])
        calculated = np.append(calculated, pair)
        radec = np.append(radec, ra[index])
        radec = np.append(radec, dec[index])
    oc = np.array([])
    for i in range(len(radec)):
        temp = radec[i] - calculated[i]
        oc = np.append(oc, temp)
  
    # deriv lists of ra and dec
    partials = np.empty((len(ra)*2,6))
    r = np.array([r2[0], r2[1], r2[2], r2dot[0], r2dot[1], r2dot[2]])

    for i in range(len(ra)):
        for j in range(6):
            deltaVec = np.array([float(0),float(0),float(0),float(0),float(0),float(0)])
            deltaVec[j] = delta # add deltaVec to r to change one term in r

            # modify r
            tempplus = np.array(r+deltaVec)
            tempr2plus = tempplus[:3]
            tempr2dotplus = tempplus[3:]
            tempminus = np.array(r-deltaVec)
            tempr2minus = tempminus[:3]
            tempr2dotminus = tempminus[3:]

            plusElements = orbitalElements(tempr2plus, tempr2dotplus, t[i])
            minusElem = orbitalElements(tempr2minus, tempr2dotminus, t[i])
            raPlus, decPlus = generateEphem(plusElements[0], plusElements[1], plusElements[2], plusElements[3], plusElements[4], plusElements[5], plusElements[6], t[i], R[i])
            raMinus, decMinus = generateEphem(minusElem[0], minusElem[1], minusElem[2], minusElem[3], minusElem[4], minusElem[5], minusElem[6], t[i], R[i])
            raderiv = (raPlus - raMinus)/(2*delta)
            decderiv = (decPlus - decMinus)/(2*delta)
            partials[i][j] = raderiv
            partials[i+len(ra)][j] = decderiv
            
    sumsderiv = np.array([])
    ocsums = np.array([])
    add = 0
    for i in range(6):
        for j in range(len(ra)*2):
            cutPartial = partials[:,i:i+1]
            add+= np.sum(oc[j]*cutPartial[j])
        ocsums = np.append(ocsums, add)
    for i in range(6):
        vals = partials[:,i:i+1]
        sumvals = np.sum(vals)
        sumsderiv = np.append(sumsderiv, sumvals)
    jMatrix = np.empty((6,6))
    for i in range(6):
        left = sumsderiv[i]
        for j in range(6):
            right = sumsderiv[j]
            jMatrix[i][j] = left*right
    xMatrix = np.linalg.lstsq(jMatrix, ocsums, rcond=None)
    xMatrix = xMatrix[0]
    return xMatrix

# inputs are lists
def monteCarlo(ra, dec, uncertaintyRa, uncertaintyDec, t, R, rootNum):
    K = 10000 # 10,000 runs through monte carlo
    start = time.time()
    for a in range(1, K+1):
        monte = open("monteCarlo.txt", "w")
        monteRA = np.array([])
        monteDec = np.array([])
        N = 1
        for i in range(len(ra)):
            sampleRa = np.random.normal(ra[i], radians(uncertaintyRa[i]), N)
            sampleDec = np.random.normal(dec[i], radians(uncertaintyDec[i]), N)
            monteRA = np.append(monteRA, sampleRa)
            monteDec = np.append(monteDec, sampleDec)
        monteRA = np.reshape(monteRA, (len(ra), N))
        monteDec = np.reshape(monteDec, (len(dec), N))
        for i in range(len(ra)):
            for j in range(N):
                cutRA = monteRA[i:i+1, j:j+1]
                cutDec = monteDec[i:i+1, j:j+1]
                monte.write(str(t[i]))
                monte.write(" ")
                monte.write(str(cutRA[0][0]))
                monte.write(" ")
                monte.write(str(cutDec[0][0]))
                monte.write(" ")
                monte.write(str(R[i][0]))
                monte.write(" ")
                monte.write(str(R[i][1]))
                monte.write(" ")
                monte.write(str(R[i][2]))
                monte.write("\n")
        print("count :", a)
        monte.close()
        runMonte("monteCarlo.txt", rootNum)
    print("Duration:  ", time.time() - start, "\n")

