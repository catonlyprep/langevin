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

import ctypes, os
from math import pi, sqrt
import numpy as np
from numpy.ctypeslib import ndpointer
import h5py

class Langevin(object):

    """       

    Parameters
    ----------

    Attributes
    ----------

    """
    def __init__(self, potential, xcell=(0,1), ycell=(0,1), initpos=(0,0), length=1, dt=1, temp=300, mass=1, gamma=1, seed=0):
        # Physical constants (2010 CODATA)
        self._e   = 1.602176565e-19  # Elementary charge
        self._ang = 1.0e-10          # Angstrom
        self._k   = 1.3806488e-23    # Boltzmann constant, J/K 
        self._amu = 1.660538921e-27  # Atomic mass unit, kg 
        
        # Time constants
        self._s   = length * sqrt(self._e / self._amu) / self._ang
        self._fs  = 1e-15 * self._s
        
        # Boltzmann's constant, eV/K
        self._kb  = self._k / self._e
        
        # Simulation parameters
        self.initpos = initpos
        self.xcell   = xcell
        self.ycell   = ycell
        self.dt      = dt * self._fs
        self.temp    = temp * self._kb
        self.mass    = mass
        self.gamma   = gamma
        self.seed    = int(seed)
        self.pot     = potential
        
        # Setup ctypes function
        libpath = os.path.dirname(os.path.abspath(__file__)) + '/liblangevin.so'
        self._Simulation = ctypes.cdll.LoadLibrary(libpath).LangevinSimulator
        self._Simulation.restype = ctypes.c_int
        self._Simulation.argtypes = [
                           ctypes.c_char_p,
                           ctypes.c_char_p,
                           ndpointer(ctypes.c_double),
                           ndpointer(ctypes.c_double),  
                           ctypes.c_double,
                           ctypes.c_double,
                           ctypes.c_double,
                           ctypes.c_double,
                           ctypes.c_int,
                           ctypes.c_int,
                           ctypes.c_int]

  
    def run(self, fname, nsteps, nsnapshots):
        """Run the Langevin simulation
        
        Parameters
        ----------
        fname : string
            Output filename for trajectory
            
        nsteps : integer
            Number of steps to simulate
            
        nsnapshots : integer
            Size between steps to return
        
        Returns
        -------
        result : integer
            0 for successful simulation
            1 for errors (needs implementing!)
        
        """
        pos = np.asarray(self.initpos)
        cell = np.asarray(self.xcell + self.ycell)                       
        result = self._Simulation(
                        fname, 
                        self.pot, 
                        cell, 
                        pos, 
                        self.dt, 
                        self.temp, 
                        self.gamma, 
                        self.mass, 
                        int(nsteps),
                        int(nsnapshots), 
                        self.seed)
        return result      

    def load(self, fname):
        """Load the output trajectory
        
        Parameters
        ----------
        fname : string
            Output filename for trajectory
        
        Returns
        -------
        trajectory : array [2, nsteps]
            The trajectory as a 2xnstep array
        
        """
        f = h5py.File(fname, 'r')
        trajectory = f['dataset'][:]
        return trajectory



