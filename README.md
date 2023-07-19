# BIMBAMBUM v1.0

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
nearby boundary. The associated code is available in the `Bubble_Dynamics` folder. The second stage computes the velocity and pressure 
fields associated with the bubble dynamics at any selected time point in the bubble lifetime. 
The associated code is available in the `Flow_Field_Quantities` folder. Both stages of the solver must be built and executed separately.
The first stage of the solver may be used independently while the second stage needs
the simualtion results of the first stage as inputs.

    BIMBAMBUM v1.0
    ├── Bubble_Dynamics         # Bubble and fluid-fluid interface dynamics
    │   ├── src                 # c++ source files
    │   ├── include             # c++ header files 
    │   ├── solver              # main function for solver execution, called at program startup.
    │   ├── config              # example config files
    │   ├── post_processing     # python scripts for rudimentary post-processing of the simulation
    │   ├── doc                 # doygen folder (compilation optional)
    │   └── CMakeLists.txt 
    ├── Flow_Field_Quantities   # Velocity and pressure fields
    │   ├── CPP_Code            # c++ source and header files
    │   ├── Python_Code         # python scripts 
    │   ├── Examples            # examples input files
    │   └── CMakeLists.txt 
    ├── license.md
    └── README.md
    
## Dependencies

The solver relies on the following dependencies.

* c++:
    * a C++ compiller with OpenMP support (tested: Clang 11.0.3 on macOS and )
    * a CMake build system (https://cmake.org/, version 3.14 or higher)
    * the Armadillo library for linear algebra (https://arma.sourceforge.net/docs.html; version tested 12.4.0).
    * the GNU scientific library for numerical computing (https://www.gnu.org/software/gsl/; version tested 2.7).
    * the BOOST library for input parsing (https://www.boost.org/; version tested 1.81.0)

* Python:
    * Python (version tested 3.11.3)
    * Numpy (version tested 1.24.3)
    * Matplotlib (version tested 3.7.1)
    
* Python and c++ binding:
    * pybind11 (https://github.com/pybind/pybind11; version tested 2.10.4)


## Building the software and running the examples

### Bubble dynamics

This portion of the solver computes the temporal evolution of the bubble and fluid-fluid interface boundaries. It may be 
run as stand-alone and can be built from the terminal with the following command lines:

```
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

The executable is installed in a `bin` folder alongside a set of `.json` configuration files provided as examples to run
the simulations. These configurations files contain the input data defining the nature of the simulation. The set of 
input data is provided in the table below:

| Parameters         | Type    | Short description                                                                  |
| -------------------| ------- | ---------------------------------------------------------------------------------- |
| Nb                 | int     | The number of elements for the bubble discretization                               |
| Ns                 | int     | The number of elements for the fluid-fluid interface discretization              |
| bubble_dynamics    | string  | The physics governing the bubble (`Rayleigh_Bubble` or `Rayleigh_Plesset_Bubble`)  |
| zeta               | double  | Buoyancy parameter                                                                 |
| gamma              | double  | Stand-off distance                                                                 |
| alpha              | double  | Density ratio of the fluid-fluid interface                                         |
| epsilon            | double  | Strength parameter                                                                 |
| k                  | double  | Specific heat ratio                                                                |
| surface_elasticity | boolean | Whether or not the interface has elasticity (`true` or `false`)                    |
| sigma_s            | double  | Fluid-fluid interface surface tension                                              |
| boundary           | string  | Nature of the initial geometry of the boundary (`from_code` must be used in v1.0)  |
| temporal_solver    | string  | Order of the temporal discretization (`RK1` or `RK2`)                              |
| dumper_filename    | string  | Name of the dumper file (stores the solution results)                              |
| delta_phi          | double  | Constant for adaptive time stepping                                                |
| filtering_freq     | int     | Filtering frequency (i.e. filtering every nth time steps)                          |
| n_threads          | int     | Number of threads used for parallel processing                                     |


The solver may then be run with the following command lines:

```
$ cd bin/
$ ./main <name of config file>.json
```

The results of the simulation are written in a `.txt` file whose name is defined by `dumper_filename`. The 10 first lines 
of this file contain information regarding key aspects of the simulation. For lines 11 and below, the first column indicates 
the time step and the second column the simulation time. Information regarding the other columns are documented in 
 `Bubble_Dynamics/include/BIM_solver.hpp` (in the method `write_solution`) and their access and use are succinctly exemplified 
in the scripts available in `Bubble_Dynamics/post_processing`.

### Flow field quantities

This portion of the solver computes the flow field quantities (velocity and pressure) associated with the 
bubble dynamics at any selected time point in the bubble lifetime. To build this portion of the code, the 
user must first downlaod pybind11 (https://github.com/pybind/pybind11) . This can be done as follow 
from the terminal with the following command lines::

```
$ cd Flow_Field_Quantities/
$ git clone https://github.com/pybind/pybind11.git
```
Once pybind11 is downloaded, the code may be built with the following command lines:
```
$ mkdir build
$ cd build
$ cmake ../
$ make
```

All files needed to run this portion of the solver are copied in the `build` folder alongside sample 
results of the first processing stage which are needed as inputs for the second processing stage. The remainder
of the inputs needed for the simulation may be changed by the user in the `inputs.py` script. The set of 
input data is provided in the table below: 

| Parameters         | Type         | Short description                                                                               |
| -------------------| ------------ | ----------------------------------------------------------------------------------------------- |
| file_name          | string       | Input file name (i.e. results files from first processing stage)                                |
| dumper_name        | string       | Output file name (no need to provide an extension)                                              |
| Nb                 | int          | The number of elements for the bubble discretization                                            |
| Ns                 | int          | The number of elements for the fluid-fluid interface discretization                             |
| zeta               | double       | Buoyancy parameter                                                                              |
| time_step          | list of int  | List of the time steps (from the first processing phase) for which the flow fields are computed |
| h_grid             | double       | Grid spacing                                                                                    |
| s_grid             | double       | Offsetting constant for shifting the domain boundaries away from the original boundaries        |
| r_max              | double       | Extend of the computational domain in the r-direction                                           |
| z_min              | double       | Extend of the computational domain in the z-direction (min) (N.B. the bubble is located at z<0) |
| z_max              | double       | Extend of the computational domain in the z-direction (max)                                     |


An example may then be run with the following command line:  

```
$ python3 main.py
```

The results of the simulation are written in a `<dumper_name><time_step>.txt` file.
This file may be read in the open source post-processing visualization engine `Paraview` (https://www.paraview.org/) with the
following procedure: 
1) open the file `<time_step>.txt` in `Paraview`, 
2) apply the `Table To Points` filter and assign the appropriate columns to the `X`, `Y` and `Z` coordinates and, 
3) apply the `Delaunay 2D` filter.

## Funding

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**BIMBAMBUM**  is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
**BIMBAMBUM**  is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details http://www.gnu.org/licenses.

