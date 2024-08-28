#ifndef UPDATE_MACRO_HPP
#define UPDATE_MACRO_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
#include <functional>
// ======================================

// =========== EIGEN INCLUDES ===========
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "../../utils/loggers/Logger.hpp"
#include "../../utils/aliases.hpp" 
#include "../PolicyTypes.hpp"
// =======================================

void update_macro(const Tensor<double, 3> &populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u);

#endif