#ifndef LLALBM_HPP
#define LLALBM_HPP

// Main Lattice objects
#include "../src/core/lattices/Lattice.hpp"
#include "../src/core/boundaries/BounceBackPolicy.hpp"
#include "../src/core/boundaries/ZouHePolicy.hpp"
#include "../src/core/collisions/BGKCollision.hpp"
#include "../src/core/initializers/InletOutletInitializer.hpp"
#include "../src/core/lattices/LatticeConfiguration.hpp"

// Utils - Logger utility
#include "../src/utils/loggers/Logger.hpp"
#include "../src/utils/readers/LatticeReader.hpp"
#include "../src/utils/aliases.hpp"
#include "../src/utils/generation/ConstructionInfo.hpp"
#include "../src/utils/generation/ConstructionTools.hpp"
#include "../src/utils/MultiDimensionalLoop.hpp"

#endif // LLALBM_HPP