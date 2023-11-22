#ifndef HH_LATTICE_MATRIX
#define HH_LATTICE_MATRIX

#include <iostream>
#include <string>
#include <vector>
#include "latticeNode.hpp"

template <std::size_t dim, std::size_t numVelocities>
class Lattice
{
private:
    std::vector<std::vector<LatticeNode<dim, numVelocities>> LBMMatrix;
public:
    Lattice(std::string matrixTemplate);
    
};






#endif