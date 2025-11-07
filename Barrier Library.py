# -*- coding: utf-8 -*-
"""-------------------------------------------------------------
Barrier Library
-------------------------------------------------------------"""

#Triangle
V = 0*x
for i in range(len(V)):
    if x[i] > -0.5 and x[i] < 0.25: #linear ramp up; peak at 0.25
        
        V[i] = V0 * (x[i] - (-0.5)) / (0.25 - (-0.5))
    elif x[i] >= 0.25 and x[i] < 1.0: #linear ramp down; peak at o.25
        
        V[i] = V0 * (1.0 - x[i]) / (1.0 - 0.25)

#Gaussian
V = 0*x
centre = 0.25 #centre position; change this to change the position
width = 0.3 #width
for i in range(len(V)):

    V[i] = V0 * np.exp(-(x[i] - centre)**2 / (2*width**2))

#Rectangle 
for i in range(len(V)):
    if x[i]>0 and x[i]<1: #width of well
        V[i]=V0

#Double Rectangle
V = 0*x
for i in range(len(V)):
    if (x[i] > -0.8 and x[i] < -0.3) or (x[i] > 0.3 and x[i] < 0.8): # (first barrier) or (second barrier). 
        V[i] = V0
