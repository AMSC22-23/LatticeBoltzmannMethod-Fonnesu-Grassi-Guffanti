cmake_minimum_required(VERSION 3.10)
project(lattice_boltzmann)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -Wno-ignored-qualifiers")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(SOURCES
    src/main.cpp
    src/lbm/boundary.cpp
    src/lbm/collisionModel.cpp
    src/lbm/lattice.cpp
    src/lbm/lattice2D.cpp
    src/lbm/latticeReader.cpp
    src/lbm/lbm.cpp
    src/lbm/velocitySet.cpp
)

set(HEADERS
    include/boundary.hpp
    include/collisionModel.hpp
    include/lattice.hpp
    include/lattice2D.hpp
    include/latticeNode.hpp
    include/latticeReader.hpp
    include/lbm.hpp
    include/velocitySet.hpp
    include/inletInitializer.hpp
)


add_executable("${PROJECT_NAME}-seq" ${SOURCES} ${HEADERS})

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
set(OPENMP_FLAGS "-fopenmp")
endif()

add_executable("${PROJECT_NAME}-omp" ${SOURCES} ${HEADERS})
target_compile_options("${PROJECT_NAME}-omp" PRIVATE ${OPENMP_FLAGS})

target_include_directories("${PROJECT_NAME}-seq" PRIVATE 
${EIGEN3_INCLUDE_DIRS}
include/
)

target_include_directories("${PROJECT_NAME}-omp" PRIVATE 
${EIGEN3_INCLUDE_DIRS}
include/
)

target_link_libraries ("${PROJECT_NAME}-seq" Eigen3::Eigen)
target_link_libraries ("${PROJECT_NAME}-omp" Eigen3::Eigen OpenMP::OpenMP_CXX)

#target_link_libraries (${PROJECT_NAME} Eigen3::Eigen)
