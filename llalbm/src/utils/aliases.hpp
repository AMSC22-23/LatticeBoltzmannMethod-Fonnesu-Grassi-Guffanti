#ifndef LLALBM_ALIASES_HPP
#define LLALBM_ALIASES_HPP

// =========== STL INCLUDES =============
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


// =======================================================================================================
//              COORDINATES, NODES and BOUNDARIES
// =======================================================================================================

/// Generic coordinate of a point
using Coord = Eigen::Index;

/// @brief Different types of boundaries, identified when the lattice is created/read from file.
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
enum InputNodeType {
    FLUID = 0,
    SOLID = 1, 
    BOUNDARY = 2, 
    INLET = 3, 
    OUTLET = 4, 
    OBSTACLE = 5
};

/// @brief Struct describing the coordinates of a node
/// @tparam d umber of dimensions, directly indicating the total number of coordinates 
template<std::size_t d>
struct Point 
{
    Eigen::Index coords[d];
};

/// @brief Struct describing a point at a boundary: in this case, it may be necessary to
/// know what type of node (wall, edge, angle) a boundary is
/// @tparam d Number of dimensions, directly indicating the total number of coordinates
template<std::size_t d>
struct BoundaryPoint
{
    Eigen::Index coords[d];
    btypes type;
};


// =======================================================================================================
// =======================================================================================================
// =======================================================================================================


// =======================================================================================================
//                  VELOCITY SETS (for 2D and 3D)
// =======================================================================================================

/// @brief D2Q9 velocity set
extern Matrix<double,9,3> D2Q9;

/// @brief D3Q19 velocity set
extern Matrix<double,19,4> D3Q19;

// =======================================================================================================
// =======================================================================================================
// =======================================================================================================


#endif