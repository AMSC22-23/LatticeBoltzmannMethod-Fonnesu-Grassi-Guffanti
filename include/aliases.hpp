#ifndef HH_LBM_ALIASES
#define HH_LBM_ALIASES

#include <vector>
#include <tuple>
#include <memory>

#include "latticeNode.hpp"


using LatticeGrid2D = std::vector<std::vector<LatticeNode<2>>>;
using BoundaryList2D = std::vector<std::tuple<LatticeNode<2>*, std::size_t, std::size_t>>;

#endif