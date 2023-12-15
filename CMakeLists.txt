cmake_minimum_required(VERSION 3.10)
project(lattice_boltzmann)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -Wno-ignored-qualifiers")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

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
)


add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

target_include_directories(${PROJECT_NAME} PRIVATE 
    ${EIGEN3_INCLUDE_DIRS}
    include/
)

target_link_libraries (${PROJECT_NAME} Eigen3::Eigen OpenMP::OpenMP_CXX)
