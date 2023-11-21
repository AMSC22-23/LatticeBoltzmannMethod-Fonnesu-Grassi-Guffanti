
#ifndef HH_LATTICE_MATRIX
#define HH_LATTICE_MATRIX

#include <iostream>
#include <string>
#include <vector>
#include "latticeNode.hpp"

class Lattice
{
private:
    std::vector<std::vector<LatticeNode>> LBMMatrix;
public:
    Lattice(std::string matrixTemplate);
    
};






#endif