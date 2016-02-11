# Langevin

**Simple Langevin Dynamics simulator for a particle in a 2D potential**

This is a simple Langevin Dynamics simulator for a particle moving in a 2D potential, for example 
to model a particle diffusing on a surface.

---

## Contents

+ [Installation](#installation)
+ [Usage](#usage)

## Installation

**Dependencies**

The program uses interpolation functions from version 2.1 of the  **[GNU Scientific Library](https://www.gnu.org/software/gsl/)**,
so please make sure you use the latest version.

**Building from source**

To build the library, unpack the source and `cd` into the unpacked directory, then type `make`:

```bash
$ tar -xzf langevin.tar.gz
$ cd langevin
$ make
```

This will generate a C++ library called `liblangevin.so`, as well as a standalone program called `langevin`.

## Usage - Python

```python
from pylangevin import Langevin

simulation = pylangevin.Langevin(options)
simulation.run(outputfilename, nsteps, nsnapshots)
```

## Usage - Standalone

```bash
$ ./langevin --args
```

## Tasks

- [x] Convert to object-oriented model
- [x] Python wrapper
- [x] Remove dependency on Armadillo
- [ ] Standalone program needs to take arguments
- [x] Setup conversion for reduced units
- [ ] Write simple test function
- [ ] Basic documentation
- [ ] Examples


