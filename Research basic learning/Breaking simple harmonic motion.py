from visual import *
from visual.graph import *

## Physical constants.
mass    = 1.0 # mass of pendulum
length  = 0.25*9.8 # lendth of pendulum
damping = 0.1 # damping coefficient
grav    = 9.8 # magnitude of graviational field

gdisplay(x=200,y=400)
fdvstime = gcurve(color=color.white)

## Differential equation for second derivative of theta.
def thetadotdot(theta, thetadot, time):

    pi = 3.14159265359

##    fd = 0.0 ## Uncomment if you're taking out driving force!

# Calculate driving force.
    period = pi/2.0
##    fd = 5.0*cos(6.28*time/period)
    if (abs(time%period) < period/2.0):
        fd = 5.0
    else:
        fd = -5.0

    fdvstime.plot(pos=(time,fd))

#    if (abs(thetadot) < 0.1 and theta > 0 and time > 1.0):
##    if (abs(thetadot) < 0.1 and time > 1.0 and abs(theta) > 0.1):
##        fd = -10.0*theta/abs(theta)
##    else:
##        fd = 0.0


# The diff eq!
    thetadotdot = -grav/length*sin(theta) - damping*thetadot + fd
    return thetadotdot

## Initialize run.
time     = 0.0
tmax     = 300
dt       = 0.005
theta    = 0.5
thetadot = -0.02

## Create visuals.
gdisplay(x=0,y=0,title = 'theta vs. time')
thetavstime = gcurve(color=color.green)
gdisplay(x=500,y=400,title = 'thetadot vs. theta')
thetadotvstheta = gcurve(color=color.green)
display(x=800,y=0)
bob = sphere(pos=(length*sin(theta),-length*cos(theta),0),
             radius = length/10.0)
rod = cylinder(pos=(0,0,0),axis=bob.pos,radius=bob.radius*0.1)
rod0 = cylinder(pos=(0,0,0),axis=bob.pos,radius=bob.radius*0.1,
                color=color.red,
                opacity = 0.25) # Reference line.

while time < tmax:
    rate(100)
    thetadot = thetadot + thetadotdot(theta, thetadot, time)*dt
    theta    = theta + thetadot*dt
    time     = time + dt
    bob.pos  = (length*sin(theta),-length*cos(theta),0)
    rod.axis = bob.pos
    thetavstime.plot(pos=(time,theta))
    thetadotvstheta.plot(pos=(theta,thetadot))
