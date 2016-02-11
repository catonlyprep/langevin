/***************************************************************************

	C++ 2D Langevin Dynamics simulator

	Author:	Tom Furnival	
	Email:	tjof2@cam.ac.uk

	Copyright (C) 2015-2016 Tom Furnival
	
	Simulation units are:
        time        : femtoseconds
        mass        : atomic mass units
        length      : Angstroms
        temperature : Kelvin  

***************************************************************************/

// Option parser
#include "optionparser/ezOptionParser.hpp"

// My libraries
#include "potentials.hpp"
#include "simulation.hpp"

// This is the Python/C interface using ctypes
//	- Is C-style for simplicity
extern "C" int LangevinSimulator(
				char *trajectory,
				char *potential,
				double *cell,
				double *pos, 
				double dt, 
				double temperature, 
				double gamma, 
				double mass, 
				int nsteps,
				int snapshot,
				int rngseed ) {
	
	std::cout<<std::string(30,'-')<<std::endl;
	std::cout<<"2D Langevin Dynamics Simulator"<<std::endl<<std::endl;
	std::cout<<"Author: Tom Furnival"<<std::endl;
	std::cout<<"Email:  tjof2@cam.ac.uk"<<std::endl;
	std::cout<<std::string(30,'-')<<std::endl<<std::endl;

	std::string PotentialFileName(potential);
	std::string TrajectoryFileName(trajectory);

	// Create potential
	Potential *Pot = new Potential();
	Pot->load(PotentialFileName);
	Pot->initialize(cell);

	// Create simulation
	Simulation *Sim = new Simulation();
	
	// Initialize with parameters
	Sim->initialize(dt, temperature, gamma, mass, nsteps, snapshot, rngseed);
	
	// Set initial position
	Sim->position(pos[0], pos[1]);
	
	// Load the potential
	Sim->potential(Pot);
	
	// Run the simulation
	Sim->run();
	
	// Save trajectory
	Sim->save(TrajectoryFileName);
	
	// Free memory
	delete Sim;
	delete Pot;
	
	// Return
	std::cout<<std::string(30,'-')<<std::endl<<std::endl;
	return 0;
}

// Main program
int main(int argc, const char * argv[]) {

	ez::ezOptionParser opt;
	opt.overview = "\n2D Langevin Dynamics Simulator";
	opt.syntax = "./langevin";
	opt.example = "./langevin\n\n";
	opt.add(
		"", 							// Default.
		0, 								// Required?
		0, 								// Number of args expected.
		0, 								// Delimiter if expecting multiple args.
		"Display usage instructions.", 	// Help description.
		"-h",     						// Flag token. 
		"-help",  						// Flag token.
		"--help", 						// Flag token.
		"--usage" 						// Flag token.
	);
	opt.add("",	1, 1, 0, "Filename for output of trajectory in HDF5 format", "-o", "-output");
	opt.add("", 1, 1, 0, "Filename for input potential in HDF5 format", "-p","-potential");	
	opt.add("0,1,0,1", 1, 4, ',', "Simulation cell size: x1,x2,y1,y2", "-c", "-cell");
	opt.add("0,0", 1, 2, ',', "Particle starting position: x,y", "-s", "-start");	
	opt.add("1", 1,	1, 0, "Friction coefficient, gamma", "-g", "-gamma");
	opt.add("1", 1,	1, 0, "Mass (a.m.u)", "-m", "-mass");
	opt.add("1", 1,	1, 0, "Length (Angstroms)", "-l", "-length");
	opt.add("1", 1,	1, 0, "Temperature (Kelvin)", "-t", "-temperature");
	opt.add("1", 1,	1, 0, "Timestep (fs)", "-d", "-timestep");
	opt.add("1", 1,	1, 0, "Number of steps", "-n", "-nsteps");
	opt.add("1", 1,	1, 0, "Interval between saved steps", "-i", "-interval");
	opt.add("0", 0,	1, 0, "Random seed", "-r", "-random");
	
	// ./langevin -o tests/test_trajectory2.h5 -p tests/test_potential.h5 -c 0,1,0,1 -s 0.5,0.0 -g 10 -m 63.456 -l 2.45 -t 500 -h 0.5 -r 1 -n 10000000 -i 10000
	
	// Check for errors
	opt.parse(argc, argv);
	if (opt.isSet("-h")) {
		std::string usage;
		opt.getUsage(usage,80,ez::ezOptionParser::INTERLEAVE);
		std::cout<<usage;
		return 1;
	}
	std::vector<std::string> badOptions;
	if(!opt.gotRequired(badOptions)) {
		for(int i=0; i < (int)badOptions.size(); ++i) {
			std::cerr<<"ERROR: Missing required option "<<badOptions[i]<<std::endl;
		}
		return 1;
	}
	if(!opt.gotExpected(badOptions)) {
		for(int i=0; i < (int)badOptions.size(); ++i) {
			std::cerr<<"ERROR: Got unexpected number of arguments for option "<<badOptions[i]<<std::endl;
		}		
		return 1;
	}
	
	// Parse options
	std::string outfile, potfile;
	std::vector<double> cellsize, start;
	double gamma, length, temp, dt, mass;
	int nsteps, nsnapshots, seed; 
	
	opt.get("-output")->getString(outfile);
	opt.get("-potential")->getString(potfile);
	opt.get("-cell")->getDoubles(cellsize);
	opt.get("-start")->getDoubles(start);
	opt.get("-gamma")->getDouble(gamma);
	opt.get("-mass")->getDouble(mass);
	opt.get("-length")->getDouble(length);
	opt.get("-temperature")->getDouble(temp);
	opt.get("-timestep")->getDouble(dt);
	opt.get("-nsteps")->getInt(nsteps);
	opt.get("-interval")->getInt(nsnapshots);
	opt.get("-random")->getInt(seed);

	///////////////////////////////////////////////////////////////////
	
	std::cout<<std::string(30,'-')<<std::endl;
	std::cout<<"2D Langevin Dynamics Simulator"<<std::endl<<std::endl;
	std::cout<<"Author: Tom Furnival"<<std::endl;
	std::cout<<"Email:  tjof2@cam.ac.uk"<<std::endl;
	std::cout<<std::string(30,'-')<<std::endl<<std::endl;
	
	// Physical constants (2010 CODATA)	
	double e_charge = 1.602176565e-19;	// Elementary charge
	double ang 		= 1e-10;			// Angstrom
	double k	   	= 1.3806488e-23;	// Boltzmann constant, J/K		
	double kb  		= k / e_charge;		// Boltzmann's constant, eV/K
	double amu 		= 1.660538921e-27;	// Atomic mass unit, kg 
	double s  		= length * sqrt(e_charge / amu) / ang;
	double fs  		= 1e-15 * s;		// Femtosecond
        
	// Unit cell sizes
	double *celldims;
	celldims = &cellsize[0];
	
	// Create potential
	Potential *Pot = new Potential();
	Pot->load(potfile);
	Pot->initialize(celldims);

	// Create simulation
	Simulation *Sim = new Simulation();
	
	// Initialize with parameters
	Sim->initialize(dt * fs, temp * kb, gamma, mass, nsteps, nsnapshots, seed);
	
	// Set initial position
	Sim->position(start[0], start[1]);
	
	// Load the potential
	Sim->potential(Pot);
	
	// Run the simulation
	Sim->run();
	
	// Save trajectory
	Sim->save(outfile);
	
	// Free memory
	delete Sim;
	delete Pot;
	
	// Return
	std::cout<<std::string(30,'-')<<std::endl<<std::endl;
	return 0;
}
