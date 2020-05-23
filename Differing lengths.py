
from visual import *
from visual.graph import *

## Create graphs for Poincare maps.
#display(x=800,y=0)
poincare_scene1 = gdisplay(title='Poincare: theta_2 vs. theta_1', # Green pendulum.
                    x=0,y=0,
                    height = 1200, width = 1200,
                    #xmin=0,xmax=1.0,
                    #ymin=0,ymax=1.0,
                    xtitle='theta_1',ytitle='theta_2')
poincare_graph1 = gdots(color=color.green)
##poincare_scene2 = gdisplay(title='Poincare: theta_4 vs. theta_3', # White pendulum.
##                    x=0,y=0,
##                    height = 1200, width = 1200,
##                    #xmin=0,xmax=1.0,
##                    #ymin=0,ymax=1.0,
##                    xtitle='theta_3',ytitle='theta_4')
##poincare_graph2 = gdots(color=color.white)

## Physical constants.
mass1   = 1.0 # mass of pendulum 1
mass2   = 1.5*mass1 # mass of pendulum 2
mtot    = mass1 + mass2 # total mass
length1 = 0.25*9.8 # lendth of pendulum 1
length2 = 1.0*length1 # length of pendulum 2
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

    energy = -mass1*grav*length1*cos1-mass2*grav*(length1*cos1+length2*cos2)
    energy += 0.5*mass1*theta1dot**2*length1**2+0.5*mass2*(theta1dot**2*length1**2+theta2dot**2*length2**2+2*length1*length2*theta1dot*theta2dot*cos12)
    
    secondDerivs = [theta1dotdot, theta2dotdot, energy] ## Store results.
    
    return secondDerivs

## Initialize run.
time     = 0.0
tmax     = 400
dt       = 0.005

theta1    = 0.50*pi ## Green double pendulum.
theta1dot = 0.0
theta2    = 0.50*pi
theta2dot = 0.0

##theta3    = 0.50*pi ## White double pendulum.
##theta3dot = 0.0
##theta4    = 0.60*pi
##theta4dot = 0.0
##
disp = 4*max(length1,length2)#,length3,length4) ## Distance between pendulae.

theta1dotlast = theta1dot
theta2dotlast = theta2dot
##theta3dotlast = theta3dot
##theta4dotlast = theta4dot
theta1dotspec = 0.0 # Specified value of theta1dot.
theta2dotspec = 0.0
##theta3dotspec = 0.0
##theta4dotspec = 0.0

## Create visuals.
#gdisplay(x=0,y=0,title = 'theta1 vs. time')
#theta1vstime = gcurve(color=color.green)
#theta3vstime = gcurve(color=color.white)
#gdisplay(x=500,y=400,title = 'theta2 vs. theta1')
#phasespace1 = gdots(color=color.green)
#phasespace2 = gdots(color=color.white)
##bob1 = sphere(pos=(length1*sin(theta1),-length1*cos(theta1),0),
##             radius = length1/10.0,color=color.green)
##rod1 = cylinder(pos=(0,0,0),axis=bob1.pos,radius=bob1.radius*0.1,color=color.green)
##bob2 = sphere(pos=bob1.pos+(length2*sin(theta2),-length2*cos(theta2),0),
##             radius = length2/10.0,color=color.green,make_trail=True)
##rod2 = cylinder(pos=bob1.pos,axis=bob2.pos-bob1.pos,radius=bob2.radius*0.1,color=color.green)
##
##bob3 = sphere(pos=(length1*sin(theta3)+disp,-length1*cos(theta3),0),
##             radius = length1/10.0,color=color.white)
##rod3 = cylinder(pos=(disp,0,0),axis=bob3.pos,radius=bob3.radius*0.1,color=color.white)
##bob4 = sphere(pos=bob3.pos+(length2*sin(theta4),-length2*cos(theta4),0),
##             radius = length2/10.0,color=color.white,make_trail=True)
##rod4 = cylinder(pos=bob4.pos,axis=bob4.pos-bob3.pos,radius=bob4.radius*0.1,color=color.white,MakeTrail=True)

#rod0 = cylinder(pos=(0,0,0),axis=bob.pos,radius=bob.radius*0.1,
#                color=color.red,
#                opacity = 0.25) # Reference line.

#while time < tmax:
while True:
    rate(200000000000)

    thetadotdot = secondDerivs(theta1, theta1dot, theta2, theta2dot, time)
    theta1dot = theta1dot + thetadotdot[0]*dt
    theta1    = theta1 + theta1dot*dt
    theta2dot = theta2dot + thetadotdot[1]*dt
    theta2    = theta2 + theta2dot*dt
    time      = time + dt
##    bob1.pos  = (length1*sin(theta1),-length1*cos(theta1),0)
##    rod1.axis = bob1.pos
##    bob2.pos  = bob1.pos+(length2*sin(theta2),-length2*cos(theta2),0)
##    rod2.pos  = bob1.pos
##    rod2.axis = bob2.pos - bob1.pos

##    thetadotdot = secondDerivs(theta3, theta3dot, theta4, theta4dot, time)
##    theta3dot = theta3dot + thetadotdot[0]*dt
##    theta3    = theta3 + theta3dot*dt
##    theta4dot = theta4dot + thetadotdot[1]*dt
##    theta4    = theta4 + theta4dot*dt
    time      = time + dt
##    bob3.pos  = (length1*sin(theta3)+disp,-length1*cos(theta3),0)
##    rod3.axis = bob3.pos-vector(disp,0,0)
##    bob4.pos  = bob3.pos+(length2*sin(theta4),-length2*cos(theta4),0)
##    rod4.pos  = bob3.pos
##    rod4.axis = bob4.pos - bob3.pos
    

##    theta1vstime.plot(pos=(time,mod(theta1,pi)))
##    theta3vstime.plot(pos=(time,mod(theta3,pi)))
##    phasespace1.plot(pos=(mod(theta2,pi),mod(theta1,pi)))
##    phasespace2.plot(pos=(mod(theta4,pi),mod(theta3,pi)))

    if ((theta1dot > theta1dotspec and theta1dotlast < theta1dotspec) or (theta1dot < theta1dotspec and theta1dotlast > theta1dotspec)):
        if (theta2dot > theta2dotspec and theta2dotlast < theta2dotspec) or (theta2dot < theta2dotspec and theta2dotlast > theta2dotspec):
            poincare_graph1.plot(pos=(mod(theta1,pi),mod(theta2,pi)))
##    if ((theta3dot > theta3dotspec and theta3dotlast < theta3dotspec) or (theta3dot < theta3dotspec and theta3dotlast > theta3dotspec)):
##        if ((theta4dot > theta4dotspec and theta4dotlast < theta4dotspec) or (theta4dot < theta4dotspec and theta4dotlast > theta4dotspec)):
##            poincare_graph2.plot(pos=(mod(theta3,pi),mod(theta4,pi)))
    
    theta1dotlast = theta1dot
    theta2dotlast = theta2dot
##    theta3dotlast = theta3dot
##    theta4dotlast = theta4dot
