#ifndef LLALBM_ALIASES_HPP
#define LLALBM_ALIASES_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
// ======================================

// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

using namespace Eigen; 

/// Generic coordinate of a point
using Coord = Eigen::Index;

/// Different types of boundaries, identified when the lattice is created.
enum btypes {
    NONE = 0, 
    TOP_WALL,
    RIGHT_WALL,
    BOTTOM_WALL,
    LEFT_WALL,
    RIGHT_TOP_C,
    RIGHT_BOTTOM_C,
    LEFT_BOTTOM_C,
    LEFT_TOP_C
};

/// @brief Type of nodes
enum InputNodeType {FLUID = 0, SOLID, BOUNDARY, INLET, OUTLET, OBSTACLE};

template<std::size_t d>
struct Point 
{
    Eigen::Index coords[d];
};


template<std::size_t d>
struct BoundaryPoint
{
    Eigen::Index coords[d];
    btypes type;
};

// using boundaryPoint = std::tuple<Coord, btypes>; 

/// Velocity Sets dim +1
/// D2Q9

extern Matrix<double,9,3> D2Q9;

/// D3Q19
extern Matrix<double,19,4> D3Q19;

#endif