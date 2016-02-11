/***************************************************************************

	C++ 2D Langevin Dynamics simulator

	Author:	Tom Furnival	
	Email:	tjof2@cam.ac.uk

	Copyright (C) 2015-2016 Tom Furnival

***************************************************************************/

#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

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

// PCG RNG
#include "pcg/pcg_random.hpp"

// HDF5 library
#include "H5Cpp.h"

class Simulation {
	public:
		Simulation() {};
		
		~Simulation() {};
		
		// Initialize the simulation with parameters	
		void initialize(double dt, double temperature, double gamma, double mass, int nsteps, int snapshot, int rngseed) {
								
			// Main simulation parameters
			h = dt;				// Time step
			kT = temperature;	// Temperature
			g = gamma;			// Friction coefficient
			m = mass;			// Mass of particle
			
			// Storage and other parameters			
			n = nsteps;			// Number of steps
			s = snapshot;		// Number of steps to store
			ncur = 0;			// Start at t=0
			scur = 0;			// Start at snapshot=0
			
			// Initial momentum
			px = 0.;
			py = 0.;
			
			// Initialize Q to zero
			Q.resize(2*(n/s)+2);
					
			// Seed the generator
			RNG = SeedRNG(rngseed);
		};
			
		// Set the initial position	
		void position(double x, double y) {
			qx = x;
			qy = y;
			Q[0] = x;
			Q[1] = y;
			return;
		}
				
		// Load the potential
		void potential(Potential *potential) {
			pot = potential;
			return;
		};
		
		// Save trajectory as HDF5
		void save(std::string fname) {
			H5::H5File TrajectoryFile(fname.c_str(), H5F_ACC_TRUNC);
			hsize_t dimsf[2];
			dimsf[0] = (n/s)+1;
			dimsf[1] = 2;
			H5::DataSpace dspace(2, dimsf);
			H5::DataType dtype(H5::PredType::NATIVE_DOUBLE);
			H5::DataSet dset = TrajectoryFile.createDataSet("dataset", dtype, dspace);
			dset.write(&Q[0], H5::PredType::NATIVE_DOUBLE);
			dset.close();
			dspace.close();
			TrajectoryFile.close();
			return;
		};
		
		// Run the simulation
		void run() {
			
			// Start the timer		
			time_start = std::chrono::steady_clock::now();		
			std::cout<<"  0% complete"<<std::endl;
		
			// Start the simulation			
			while(ncur < n) {
				// Integrate					
				BAOAB();				
				ncur++;
				scur++;
				
				// Write positions to array
				if( scur % s == 0 ) {
					Q[2*(scur/s)] = qx;
					Q[2*(scur/s)+1] = qy;					
				}
							
				// Report progress
				if( ncur % (n/10) == 0 ) {
					std::cout<<"  "<<(int)(100.*(double)ncur/(double)n)<<"% complete"<<std::endl;
				}
			}
					
			// Report timing
			time_end = std::chrono::steady_clock::now();
			run_time = (std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count()/1E6);
			std::cout<<std::endl<<"Total time: "<<std::setprecision(5)<<run_time<<" seconds"<<std::endl<<std::endl;	
		
			return;		
		};		
	
	private:
		// Random numbers
		pcg64 RNG;
		std::normal_distribution<double> NormalDistribution {0, 1};
		
		// Timing
		std::chrono::time_point<std::chrono::steady_clock> time_start, time_end;
		
		// Parameters
		int n, s, ncur, scur;
		
		// Parameters
		double kT, g, h, m, run_time;
		
		// Position and momentum for export
		std::vector<double> Q;

		// BAOAB paramters
		std::pair<double, double> qh, ph, qnp1, pnp1, fn;
		
		// Potential
		Potential* pot;
		
		// Current atom position, momenta and RNG
		double qx, qy, px, py, rx, ry;
		
		// Seed random number generator
		pcg64 SeedRNG(int seed) {
			// Check for user-defined seed
			if(seed > 0) {	
				return pcg64(seed);
			}
			else {				
				// Initialize random seed
				pcg_extras::seed_seq_from<std::random_device> seed_source;
				return pcg64(seed_source);
			}    
		};
		
		// BAOAB integrator
		void BAOAB() {
			// Draw random numbers
			rx = NormalDistribution(RNG);
			ry = NormalDistribution(RNG);
					
			// Exp[-gamma * dt]
			double c1 = std::exp(-1*g*h);
			
			// Sqrt[kbT * (1-c^2) * m]
			double c3 = std::sqrt(kT * (1-c1*c1) * m);
	
			// First force calculation
			fn.first  = pot->f_dx(qx, qy);
			fn.second = pot->f_dy(qx, qy);

			// Momentum half-step "B"
			ph.first  = px + (h/2)*fn.first;
			ph.second = py + (h/2)*fn.second;

			
			// Position half-step "A"
			qh.first  = qx + (h/(2*m))*ph.first;
			qh.second = qy + (h/(2*m))*ph.second;	

			// Weak solve of Ornstein-Uhlenbeck process "O"	
			ph.first  = c1*ph.first  + c3*rx;
			ph.second = c1*ph.second + c3*ry;
	
			// Position half-step "A"
			qnp1.first  = qh.first  + (h/(2*m))*ph.first;
			qnp1.second = qh.second + (h/(2*m))*ph.second;			
	
			// New force calculation
			fn.first  = pot->f_dx(qnp1.first, qnp1.second);
			fn.second = pot->f_dy(qnp1.first, qnp1.second);
				
			// Momentum half-step "B"
			pnp1.first  = ph.first  + (h/2)*fn.first;
			pnp1.second = ph.second + (h/2)*fn.second;
			
			// Update with new positions and momenta
			qx = qnp1.first;
			qy = qnp1.second;
			px = pnp1.first;
			py = pnp1.second;

			return;
		}	
};

#endif
