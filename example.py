"""
    C++ 2D Langevin Dynamics simulator

    Author:	Tom Furnival	
    Email:	tjof2@cam.ac.uk

    Copyright (C) 2015-2016 Tom Furnival
    
    Simulation units are:
        time        : femtoseconds
        mass        : atomic mass units
        length      : Angstroms
        temperature : Kelvin    
"""

import os, sys, warnings
import matplotlib.pyplot as plt
import numpy as np
from pylangevin import Langevin

simulation = Langevin(
        potential='tests/test_potential.h5',
        xcell=(0.,1.),
        ycell=(0.,1.),
        initpos=(0.5, 0.0),
        length=2.45,
        dt=0.5,
        temp=500, 
        mass=63.456, 
        gamma=10,
        seed=1
)
simulation.run('tests/test_trajectory.h5', 1e7, 1e4)

# Visualize the trajectory
trajectory = simulation.load('tests/test_trajectory.h5')

minTx = np.floor(np.amin(trajectory[:,0]))
maxTx = np.ceil(np.amax(trajectory[:,0]))
minTy = np.floor(np.amin(trajectory[:,1]))
maxTy = np.ceil(np.amax(trajectory[:,1]))

plt.figure()
plt.plot(trajectory[:,0],trajectory[:,1])
plt.xlim([minTx,maxTx])
plt.ylim([minTy,maxTy])
plt.show()

