[![Build Status](https://travis-ci.org/tjof2/langevin.svg?branch=master)](https://travis-ci.org/tjof2/langevin)

# Langevin

**Simple Langevin Dynamics simulator for a particle in a 2D potential**

This is a simple Langevin Dynamics simulator for a particle moving in a 2D potential, for example 
to model a particle diffusing on a surface.

---

## Installation

#### Dependencies

The Langevin simulator uses 2D interpolation functions from version 2.1 of the  **[GNU Scientific Library](https://www.gnu.org/software/gsl/)**,
so please make sure you use the latest version. It also uses the [HDF5 format](https://www.hdfgroup.org/HDF5/) to load and save potentials
and trajectories.

#### Building from source

To build the library, unpack the source and `cd` into the unpacked directory, then type `make`:

```bash
$ tar -xzf langevin.tar.gz
$ cd langevin
$ make
```

This will generate a C++ library called `liblangevin.so`, as well as a standalone program called `langevin`.

## Usage 

### Python

```python
from pylangevin import Langevin

simulation = pylangevin.Langevin(options)
simulation.run(outputfilename, nsteps, nsnapshots)
```

#### Standalone

The single-threaded standalone program is useful for running on HPC systems, for example when you 
want to carry out repeated runs with the same set of parameters. The program takes the following
required command-line arguments.

```bash
$ ./langevin -o tests/test_trajectory2.h5 -p tests/test_potential.h5 \
             -c 0.,1.,0.,1. -s 0.5,0.0 -g 10 -m 63.456 -l 2.45 \
             -t 500 -d 0.5 -n 10000000 -i 10000
```

You can also specify a random seed using the flag `-r [seed]`.

#### Generating your own potentials

Will be expanded on here.

## Tasks

- [ ] Write simple test function
- [ ] Basic documentation
- [ ] Examples


