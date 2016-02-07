from __future__ import division
from numpy import *
from scipy.integrate import odeint

import pylab

dt = 0.01

def eom_odeint(state, t):
	r1   = array([state[0], state[1]])
	r2   = array([state[4], state[5]])
	r21  = r2 - r1
	
	v1   = array([state[2], state[3]])
	v2   = array([state[6], state[7]])
	
	hr1  = hypot(r1[0], r1[1])
	hr2  = hypot(r2[0], r2[1])
	hr21 = hypot(r21[0], r21[1])
	
	a1 = -2*r1/(hr1**3) - r21/(hr21**3)
	a2 = -2*r2/(hr2**3) + r21/(hr21**3)
	
	return array([v1[0], v1[1], a1[0], a1[1], v2[0], v2[1], a2[0], a2[1]])
	
def angular_momentum_odeint(state):
	x1, y1, vx1, vy1, x2, y2, vx2, vy2 = state
	return (x1*vy1-y1*vx1) + (x2*vy2-y2*vx2)
	
def potential_energy(state):
	x1, y1, vx1, vy1, x2, y2, vx2, vy2 = state
	r1   = array([state[0], state[1]])
	r2   = array([state[4], state[5]])
	r21  = r2 - r1
	hr1  = hypot(r1[0], r1[1])
	hr2  = hypot(r2[0], r2[1])
	hr21 = hypot(r21[0], r21[1])
	a1 = -2*r1/(hr1**3) + r21/(hr21**3)
	a2 = -2*r2/(hr2**3) - r21/(hr21**3)
	P1 = - a1 * (vx1**2 + vy1**2)**0.5 * dt
	P2 = - a2 * (vx2**2 + vy2**2)**0.5 * dt
	return P1, P2

def energy_odeint(y):
	r1   = array([y[0], y[1]])
	v1   = array([y[2], y[3]])
	r2   = array([y[4], y[5]])
	v2   = array([y[6], y[7]])
	r21  = r2 - r1
	hr1  = hypot(r1[0], r1[1])
	hv1  = hypot(v1[0], v1[1])
	hr2  = hypot(r2[0], r2[1])
	hv2  = hypot(v2[0], v2[1])
	hr21 = hypot(r21[0], r21[1])
	
	P1, P2 = potential_energy(y)
	
	return 0.5*hv1**2 + 0.5*hv2**2 + (P1[0]**2 + P1[1]**2)**0.5 + (P2[0]**2 + P2[1]**2)**0.5 
	
y0 = array([1.4, 0.0, 0, 0.86, -1, 0.0, 0.0, -1]) #In form x1, y1, vx1, vy1...

timesn = linspace(0.0, 100.0, int(100.0 / dt))
resn = odeint(eom_odeint,y0,timesn)

pylab.figure(1)
pylab.plot([e[0] for e in resn], [e[1] for e in resn]) # position of 1
pylab.plot([e[4] for e in resn], [e[5] for e in resn]) # position of 2

pylab.figure(2)
pylab.plot([energy_odeint(e) for e in resn])
pylab.ylabel('Total Energy')

pylab.figure(3)
pylab.plot([angular_momentum_odeint(e) for e in resn])
pylab.ylabel('Angular Momentum')