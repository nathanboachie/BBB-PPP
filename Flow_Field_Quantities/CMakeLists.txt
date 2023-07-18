cmake_minimum_required(VERSION 3.18)

project(Flow_Field_Quantities)

add_subdirectory(pybind11-2.10.4)

pybind11_add_module(flow_potential
			CPP_Code/velocity_potential.cpp
			CPP_Code/cubic_spline.cpp)


# finding armadillo library and include path
# ------------------------------------------

find_library(ARMADILLO_LIBRARY libarmadillo.12.4.0.dylib /usr/local/Cellar/armadillo/12.4.0/lib/)
target_link_libraries (flow_potential PRIVATE ${ARMADILLO_LIBRARY})

find_path(ARMADILLO_INCLUDE armadillo /usr/local/Cellar/armadillo/12.4.0/include/)
include_directories(${ARMADILLO_INCLUDE})


# finding GSL library and include path
# -----------------------------------

find_library(GSL_LIBRARY libgsl.27.dylib /usr/local/Cellar/gsl/2.7.1/lib/)
target_link_libraries (flow_potential PRIVATE ${GSL_LIBRARY})

