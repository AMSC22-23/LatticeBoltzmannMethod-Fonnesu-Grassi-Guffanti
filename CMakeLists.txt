cmake_minimum_required(VERSION 3.10)

# Project name and version
project(LBM_PROJECT VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# find required packages
# find MPI
find_package(MPI REQUIRED)
# find OpenMP
cmake_policy(SET CMP0074 NEW)
find_package(OpenMP REQUIRED)
# find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE HINTS $ENV{mkEigenPrefix})


# Include llalbm library
add_subdirectory(llalbm)

# Include tests directory
add_subdirectory(test)