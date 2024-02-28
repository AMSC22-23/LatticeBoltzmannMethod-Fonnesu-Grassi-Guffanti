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

// =========== LLALBM INCLUDES ===========
#include "Logger.hpp"
// =======================================


using namespace Eigen; 
/// Generic d-dimensional point described by its set of coordinates
template<std::size_t d>
using Point = Matrix<std::size_t, d, 1>;

enum btypes {TOP_WALL, RIGHT_WALL, BOTTOM_WALL, LEFT_WALL, RIGHT_TOP_C, RIGHT_BOTTOM_C, LEFT_BOTTOM_C, LEFT_TOP_C };
template<std::size_t d>
using boundaryPoint = std::tuple<Point<d>, btypes>; 

/// Velocity Sets dim +1
/// D2Q9

extern Matrix<double,9,3> D2Q9;

/// D3Q19
extern Matrix<double,19,4> D3Q19;
