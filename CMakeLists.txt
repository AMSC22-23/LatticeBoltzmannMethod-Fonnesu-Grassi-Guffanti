cmake_minimum_required(VERSION 3.10)
project(lattice_boltzmann)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -Wno-ignored-qualifiers")
include_directories("include/")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
file(GLOB_RECURSE SOURCES src/*.cpp)


add_executable(${PROJECT_NAME} ${SOURCES})
