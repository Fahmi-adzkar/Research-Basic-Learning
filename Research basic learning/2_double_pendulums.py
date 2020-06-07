import numpy


## Physical constants.
mass1   = 1.0 # mass of pendulum 1
mass2   = mass1 # mass of pendulum 2
mtot    = mass1 + mass2 # total mass
length1 = 0.25*9.8 # lendth of pendulum 1
length2 = length1 # length of pendulum 2
damping = 0.1 # damping coefficient
grav    = 9.8 # magnitude of graviational field
pi      = 3.14159265359

## Differential equation for second derivative of theta.
def secondDerivs(theta1, theta1dot, theta2, theta2dot, time):

# Calculate frequently used stuff.
    cos1  = cos(theta1)
    cos2  = cos(theta2)
    sin1  = sin(theta1)
    sin2  = sin(theta2)
    sin12 = sin(theta1 - theta2)
    cos12 = cos(theta1 - theta2)

# The diff eqs!

    A = mtot*length1
    B = mass2*length2*cos12
    C = -mass2*length2*theta2dot**2*sin12-mtot*grav*sin1
    D = length1/length2*cos12
    E = (length1*theta1dot**2*sin12-grav*sin2)/length2

    theta1dotdot = (C-B*E)/(A-B*D)
    theta2dotdot = E-D*theta1dotdot

### BAD! Do not use.
##    theta1dotdot = -mass2*length1*theta1dot**2*sin12*cos12
##    theta1dotdot += grav*mass2*sin2*cos12
##    theta1dotdot += -mass2*length2*theta2dot**2*sin12
##    theta1dotdot += -mtot*grav*sin1
##    theta1dotdot = theta1dotdot/(length1*mtot-mass2*length1*cos12**2)
##
##    theta2dotdot = mass2*length2*theta2dot**2*sin12*cos12
##    theta2dotdot += grav*mtot*sin1*cos12
##    theta2dotdot += +mtot*length1*theta1dot**2*sin12
##    theta2dotdot += -mtot*grav*sin2
##    theta2dotdot = theta1dotdot/(length2*mtot-mass2*length2*cos12**2)

    energy = -mass1*grav*length1*cos1-mass2*grav*(length1*cos1+length2*cos2)
    energy += 0.5*mass1*theta1dot**2*length1**2+0.5*mass2*(theta1dot**2*length1**2+theta2dot**2*length2**2+2*length1*length2*theta1dot*theta2dot*cos12)
    
    secondDerivs = [theta1dotdot, theta2dotdot, energy] ## Store results.
    
    return secondDerivs

## Initialize run.
time     = 0.0
tmax     = 400
dt       = 0.005

theta1    = 0.50*pi ## Green double pendulum.
theta1dot = 0
theta2    = 0.50*pi
theta2dot = 0

theta3    = 0.50*pi ## White double pendulum.
theta3dot = 0
theta4    = 0.60*pi
theta4dot = 0

disp = 4*max(length1,length2)#,length3,length4) ## Distance between pendulae.

## Create visuals.
gdisplay(x=0,y=0,title = 'theta1 vs. time')
theta1vstime = gcurve(color=color.green)
theta3vstime = gcurve(color=color.white)
gdisplay(x=500,y=400,title = 'theta2 vs. theta1')
phasespace1 = gdots(color=color.green)
phasespace2 = gdots(color=color.white)
display(x=800,y=0)
bob1 = sphere(pos=(length1*sin(theta1),-length1*cos(theta1),0),
             radius = length1/10.0,color=color.green)
rod1 = cylinder(pos=(0,0,0),axis=bob1.pos,radius=bob1.radius*0.1,color=color.green)
bob2 = sphere(pos=bob1.pos+(length2*sin(theta2),-length2*cos(theta2),0),
             radius = length2/10.0,color=color.green,make_trail=True)
rod2 = cylinder(pos=bob1.pos,axis=bob2.pos-bob1.pos,radius=bob2.radius*0.1,color=color.green)

bob3 = sphere(pos=(length1*sin(theta3)+disp,-length1*cos(theta3),0),
             radius = length1/10.0,color=color.white)
rod3 = cylinder(pos=(disp,0,0),axis=bob3.pos,radius=bob3.radius*0.1,color=color.white)
bob4 = sphere(pos=bob3.pos+(length2*sin(theta4),-length2*cos(theta4),0),
             radius = length2/10.0,color=color.white,make_trail=True)
rod4 = cylinder(pos=bob4.pos,axis=bob4.pos-bob3.pos,radius=bob4.radius*0.1,color=color.white,MakeTrail=True)

#rod0 = cylinder(pos=(0,0,0),axis=bob.pos,radius=bob.radius*0.1,
#                color=color.red,
#                opacity = 0.25) # Reference line.

while time < tmax:
    rate(200)

    thetadotdot = secondDerivs(theta1, theta1dot, theta2, theta2dot, time)
    theta1dot = theta1dot + thetadotdot[0]*dt
    theta1    = theta1 + theta1dot*dt
    theta2dot = theta2dot + thetadotdot[1]*dt
    theta2    = theta2 + theta2dot*dt
    time      = time + dt
    bob1.pos  = (length1*sin(theta1),-length1*cos(theta1),0)
    rod1.axis = bob1.pos
    bob2.pos  = bob1.pos+(length2*sin(theta2),-length2*cos(theta2),0)
    rod2.pos  = bob1.pos
    rod2.axis = bob2.pos - bob1.pos

    thetadotdot = secondDerivs(theta3, theta3dot, theta4, theta4dot, time)
    theta3dot = theta3dot + thetadotdot[0]*dt
    theta3    = theta3 + theta3dot*dt
    theta4dot = theta4dot + thetadotdot[1]*dt
    theta4    = theta4 + theta4dot*dt
    time      = time + dt
    bob3.pos  = (length1*sin(theta3)+disp,-length1*cos(theta3),0)
    rod3.axis = bob3.pos-vector(disp,0,0)
    bob4.pos  = bob3.pos+(length2*sin(theta4),-length2*cos(theta4),0)
    rod4.pos  = bob3.pos
    rod4.axis = bob4.pos - bob3.pos
    

    theta1vstime.plot(pos=(time,mod(theta1,pi)))
    theta3vstime.plot(pos=(time,mod(theta3,pi)))
    phasespace1.plot(pos=(mod(theta2,pi),mod(theta1,pi)))
    phasespace2.plot(pos=(mod(theta4,pi),mod(theta3,pi)))
