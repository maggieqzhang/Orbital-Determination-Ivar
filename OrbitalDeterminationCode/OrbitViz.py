# Orbit visualization of 1929 SH (Ivar)

from vpython import*

# 6 orbital elements and a constant epsilon of Asteroid
a = 1.8036094125745041
e = 0.37924137179443357
M = radians(353.27849928937746)
Oprime = radians(133.50935269704902)
iprime = radians(8.197602589854448)
wprime = radians(167.7626725606778)
epsilon = radians(23.4352)

# guess of E using Kepler's Equation
def solvekep(M):
    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    Mtrue = M
    while abs(Mguess - Mtrue) > 1e-004:
        Mguess = Eguess - e*sin(Eguess)
        Eguess = Eguess - (Eguess-e*sin(Eguess)-Mtrue) / (1-e*cos(Eguess))
    return Eguess

# function to rotate the vector using numpy and rotate function
def rotateCartesian(vect, wprime, iprime, Oprime):
    vect = rotate(vect, angle = (-wprime), axis = vector(0,0,1))
    vect = rotate(vect, angle = (iprime), axis = vector(1,0,0))
    vect = rotate(vect, angle = (Oprime), axis = vector(0,0,1))
    vect = rotate(vect, angle = (epsilon), axis = vector(1,0,0))
    return vector(vect)

# defining constants and initial conditionsd
mu = 1
time = 0
dt = .05
period = 2*pi*sqrt(a**3/mu)
Mtrue = M                                                               
Etrue = solvekep(Mtrue)
rCartesian = vector(a*cos(Etrue) - a*e, a*sin(Etrue)*sqrt(1-e**2), 0)
r1ecliptic = rotateCartesian(rCartesian, wprime, iprime, Oprime)
asteroid = sphere(pos=r1ecliptic*150,radius=15, color=color.white)
asteroid.trail = curve(color=color.white)
asteroidLabel = label(pos = asteroid.pos + vector(40,30,30), text = "Ivar", height = 12, space = asteroid.radius, font = "sans")
sun = sphere(pos=vector(0,0,0), radius=50, color=color.yellow)
sunLabel = label(pos = vector(60,60,60), text = "Sun", space = sun.radius, height = 15, border = 6, font = "sans", color = color.yellow)
    
# 6 orbital elements and a constant epsilon of Earth
aE = 1.00000011	
eE = 0.01671022
ME = radians(0)
OprimeE = radians(-11.26064)
iprimeE = radians(0.00005)
wprimeE = radians(102.94719)

periodE = 2*pi*sqrt(aE**3/mu)
MtrueE = ME                                                               
EtrueE = solvekep(MtrueE)
rCartesianE = vector(aE*cos(EtrueE) - aE*eE, aE*sin(EtrueE)*sqrt(1-eE**2), 0)
r1eclipticE = rotateCartesian(rCartesianE, wprimeE, iprimeE, OprimeE)
earth = sphere(pos=r1eclipticE*150,radius=20, color=color.blue)
earth.trail = curve(color=color.white)
earthLabel = label(pos = earth.pos + vector(40,30,30), text = "Earth", height = 12, space = earth.radius, font = "sans")

# moving the Earth and Asteroid
while True:
    rate(60)
    time = time + dt
    Mtrue = 2*pi*(time)/period + M
    Etrue = solvekep(Mtrue)
    rCartesian = vector(a*cos(Etrue) - a*e, a*sin(Etrue)*sqrt(1-e**2), 0)
    r1ecliptic = rotateCartesian(rCartesian, wprime, iprime, Oprime)
    asteroid.pos = r1ecliptic*150
    asteroidLabel.pos = asteroid.pos + vector(40,30,30)
    asteroid.trail.append(pos=asteroid.pos)  

    MtrueE = 2*pi*(time)/periodE + ME
    EtrueE = solvekep(Mtrue)
    rCartesianE = vector(aE*cos(EtrueE) - aE*eE, aE*sin(EtrueE)*sqrt(1-eE**2), 0)
    r1eclipticE = rotateCartesian(rCartesianE, wprimeE, iprimeE, OprimeE)
    earth.pos = r1eclipticE*150
    earthLabel.pos = earth.pos + vector(40,30,30)
    earth.trail.append(pos=earth.pos)  
