/***************************************************************************

	C++ 2D Langevin Dynamics simulator

	Author:	Tom Furnival
	Email:	tjof2@cam.ac.uk

	Copyright (C) 2015-2019 Tom Furnival

	Simulation units are:
        time        : femtoseconds
        mass        : atomic mass units
        length      : Angstroms
        temperature : Kelvin

	This file is part of langevin.

	langevin is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	langevin is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with langevin.  If not, see <http://www.gnu.org/licenses/>.

***************************************************************************/

#ifndef _POTENTIALS_HPP_
#define _POTENTIALS_HPP_

// C++ headers
#include <algorithm>
#include <chrono>
#include <ctime>
#include <random>
#include <string>
#include <utility>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>

// GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

// HDF5 library
#include "H5Cpp.h"

// Interpolated 2D potential
class Potential {
public:
  Potential(){};

  ~Potential() {
    gsl_spline2d_free(spline);
    gsl_interp_accel_free(xacc);
    gsl_interp_accel_free(yacc);
    free(za);
  };

  // Load HDF5 potential
  void load(std::string fname) {
    H5::H5File PotentialFile(fname.c_str(), H5F_ACC_RDONLY);
    H5::DataSet dset = PotentialFile.openDataSet("/dataset");
    H5::DataSpace dspace = dset.getSpace();
    hsize_t rank;
    hsize_t dims[2];
    rank = dspace.getSimpleExtentDims(dims, NULL);
    if (rank != 2) {
      std::cout << "Potential should be 2D!" << std::endl;
    }
    nx = dims[0];
    ny = dims[1];
    potential = (double *)malloc(nx * ny * sizeof(double));

    dset.read(potential, H5::PredType::NATIVE_DOUBLE);
    dset.close();
    dspace.close();
    PotentialFile.close();
  }

  // Initializes the grid with the interpolation function
  void initialize(double *cell) {
    xmin = cell[0];
    xmax = cell[1];
    ymin = cell[2];
    ymax = cell[3];

    // Get cell size
    xscale = xmax - xmin;
    yscale = ymax - ymin;

    // Setup grid point arrays in range	of cell
    double xa[nx], ya[ny];
    for (int i = 0; i < nx; i++) {
      xa[i] = xmin + xscale * (double)i / (nx - 1);
    }
    for (int i = 0; i < ny; i++) {
      ya[i] = ymin + yscale * (double)i / (ny - 1);
    }

    // Set final array
    za = (double *)malloc(nx * ny * sizeof(double));

    // Setup GSL spline function
    spline = gsl_spline2d_alloc(T, nx, ny);
    xacc = gsl_interp_accel_alloc();
    yacc = gsl_interp_accel_alloc();

    // Loop over loaded grid, make negative now
    // for force calculations
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        gsl_spline2d_set(spline, za, i, j, -1 * potential[j * ny + i]);
      }
    }

    // Initialize spline function
    gsl_spline2d_init(spline, xa, ya, za, nx, ny);

    return;
  };

  // For each of these functions, use periodic boundary conditions
  // Potential function
  double f(double x, double y) {
    x1 = PeriodicCell(x, xmin, xmax, xscale);
    y1 = PeriodicCell(y, ymin, ymax, yscale);
    return gsl_spline2d_eval(spline, x1, y1, xacc, yacc);
  };

  // df/dx of potential
  double f_dx(double x, double y) {
    x1 = PeriodicCell(x, xmin, xmax, xscale);
    y1 = PeriodicCell(y, ymin, ymax, yscale);
    // Debugging PBCs
    // std::cout<<"("<<x<<","<<y<<"),("<<x1<<","<<y1<<")"<<std::endl;
    return gsl_spline2d_eval_deriv_x(spline, x1, y1, xacc, yacc);
  };

  // df/dy of potential
  double f_dy(double x, double y) {
    x1 = PeriodicCell(x, xmin, xmax, xscale);
    y1 = PeriodicCell(y, ymin, ymax, yscale);
    return gsl_spline2d_eval_deriv_y(spline, x1, y1, xacc, yacc);
  };

private:
  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  double *za, *potential;
  int nx, ny;
  double x1, y1;
  double xmin, xmax, ymin, ymax, xscale, yscale;
  gsl_spline2d *spline;
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;

  inline double PeriodicCell(double val, double min, double max, double scale) {
    if (val < min) {
      return val + scale * std::floor((max - val) / scale);
    } else if (val > max) {
      return val - scale * std::floor((val - min) / scale);
    } else {
      return val;
    }
  }
};

#endif
