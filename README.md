# BIMBAMBUM 1.0

## BIMBAMBUM

**BIMBAMBUM** stands for Boundary Integral Method for Bubble Analysis and Modeling 
in Bounded and Unbounded Media. It is a potential flow solver mainly written in C++ intended at modeling the behaviour 
 of a single cavitation bubble in the vicinity of an initially flat fluid-fluid interface. The software solves the
 Laplace equation using a boundary integral formulation in axisymmetric coordinates. The code has the vocation to allow 
 research teams, students and any other interested person to quickly get insights into the behaviour of cavitation bubbles
and hopefully to provide a basis for an accelerated and collaborative developments of this solver.


## What can I do with the solver?

With **BIMBAMBUM** you can solve:

* The dynamics of a cavitation bubble near a fluid-fluid interface of any density ratio.
* The pressure and velocity fields associated with the bubble dynamics.

## Repository Structure

BIMBAMBUM is a two-stage numerical flow solver. The first stage computes the time evolution of the surfaces of the bubble and 
nearby boundary. The associated code is available in the 'Bubble_Dynamics' folder. The second stage computes the velocity and pressure 
fields associated with the bubble dynamics at any selected time point in the bubble lifetime. 
The associated code is available in the 'Flow_Field_Quantities' folder. Both stages of the solver must be executed, built and compiled separately.
The first stage of the solver may be used independantly while the second stage needs
the simualtion results of the first stage as inputs.

    BIMBAMBUM
    ├── Bubble_Dynamics         # Bubble and fluid-fluid interface dynamics
    │   ├── src                 # c++ source files
    │   ├── include             # c++ header files 
    │   ├── solver              # main function for solver execution, called at program startup.
    │   ├── config              # example config files
    │   ├── doc                 # doygen folder (compilation optional)
    │   └── CMakeLists.txt 
    ├── Flow_Field_Quantities   # Velocity and pressure fields
    │   ├── CPP_Code            # c++ source and header files
    │   ├── Python_Code         # python scripts 
    │   └── CMakeLists.txt 
    ├── license.md
    └── README.md
    
## Dependencies

The solver relies on the following dependencies.

* c++:
    * a C++ compiller (tested: Clang 11.0.3 on macOS and )
    * a CMake build system (https://cmake.org/, version 3.14 or higher)
    * the Armadillo library for linear algebra (https://arma.sourceforge.net/docs.html; version tested 12.4.0).
    * the GNU scientific library for numerical computing (https://www.gnu.org/software/gsl/; version tested 2.7).
    * the OpenMP library for multiprocessing (https://www.openmp.org/).
    * the BOOST library for input parsing (https://www.boost.org/; version tested 1.81.0)

* Python:
    * Python (version tested 3.11.3)
    * Numpy (version tested 1.24.3)
    * Matplotlib (version tested 3.7.1)
    
* Python and c++ binding:
    * pybind11 (https://github.com/pybind/pybind11; version tested 2.10.4)


## Building the code
### Bubble dynamics

This portion of the solver computes the temporal evolution of the bubble and fluid-fluid interafce bounadries. It may be 
run as stand-alone and is built as follow:
```sh
$ cd Bubble_Dynamics/
$ mkdir build
$ cd build
$ cmake -DBUILD_OPENMP=<OpenMP option> -DBUILD_DOXYGEN=<Doxygen option> ../
$ make install       # optional -> install the solver as a library
$ make solver
$ make doxydoc       # optional -> write the doxygen documentation
```
where `OpenMP option` and `Doxygen option` are either 0 (false) or 1 (true) to enable parallel processing
and the built of the Doxygen documentation, respectively.

### Flow field quantities

This portion of the solver computes the flow field quantities (velocity and pressure) associated with the 
bubble dynamics at any selected time point in the bubble lifetime. To run this portion of the code, the 
user must first downlaod pybind11 (https://github.com/pybind/pybind11) . This can be done as follow 
from:

```
$ cd Flow_Field_Quantities
$ git clone https://github.com/pybind/pybind11.git
```

This portion of the code is writtne in Python, yet the performance critical tasks are written in C++ and must
be compilled. This can be done with the following command sequence conducted in the
working directory:
```
$ mkdir build
$ cd build
$ cmake ../
$ make
```

```sh
$ python main.py
```

## Running the examples
The example (such as `main.cc`) file can be built with
```sh
$ cd build
$ make examples
```

## Sample results

### Bubble dynamics
### Flow field quantities

## Funding

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**BIMBAMBUM**  is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
**BIMBAMBUM**  is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details http://www.gnu.org/licenses.

