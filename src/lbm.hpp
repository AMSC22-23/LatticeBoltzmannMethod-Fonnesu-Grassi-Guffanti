#ifndef HH_LBM
#define HH_LBM

#include <iostream>
#include <string>
#include <vector>
#include "lattice.hpp"
#include "collisionModel.hpp"
#include "velocitySet.hpp"

class lbm
{
private:
    
    VelocitySet velocitySet;

    // TODO: aggiungere tutte le costanti del caso
    int const tao;
    int const deltaT;



    /*
        function called by constructor setting all inital state
    */
    void loadInitalState();

public:
    /*
        Costructor taking input data from file
    */
    lbm(std::string input_file);

    /*
        Solving function, taked in input the iteration count t, creates files for every state of the matrix T and returns 0 if some error occured, 1 otherwise
    */
    int compute(std::size_t t);

};

#endif