# -*- coding: utf-8 -*-
# Copyright 2015-2019 Tom Furnival
#
# This file is part of langevin.
#
# langevin is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# langevin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with langevin.  If not, see <http://www.gnu.org/licenses/>.

import os, sys, warnings
import matplotlib.pyplot as plt
import numpy as np
from pylangevin import Langevin

if __name__ == "__main__":
    simulation = Langevin(
        potential="tests/test_potential.h5",
        xcell=(0.0, 1.0),
        ycell=(0.0, 1.0),
        initpos=(0.5, 0.0),
        length=2.45,
        dt=0.5,
        temp=500,
        mass=63.456,
        gamma=10,
        seed=1,
    )
    simulation.run("tests/test_trajectory.h5", 1e7, 1e4)

    # Visualize the trajectory
    trajectory = simulation.load("tests/test_trajectory.h5")

    minTx = np.floor(np.amin(trajectory[:, 0]))
    maxTx = np.ceil(np.amax(trajectory[:, 0]))
    minTy = np.floor(np.amin(trajectory[:, 1]))
    maxTy = np.ceil(np.amax(trajectory[:, 1]))

    plt.figure()
    plt.plot(trajectory[:, 0], trajectory[:, 1])
    plt.xlim([minTx, maxTx])
    plt.ylim([minTy, maxTy])
    plt.show()
