#ifndef HH_LBM
#define HH_LBM

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "lattice.hpp"
#include "lattice2D.hpp"
#include "collisionModel.hpp"
#include "velocitySet.hpp"

class lbm
{
private:
    
    VelocitySet velocitySet;

    std::unique_ptr<Lattice> lattice_ptr;
    // TODO: aggiungere tutte le costanti del caso
    double const tao;
    double const deltaT;



    /*
        function called by constructor setting all inital state
    */
    void load_initial_state();

public:
    /*
        Costructor taking input data from file
    */
    lbm(size_t D, int Q, const std::string& input_file_path);

    /*
        Solving function, taked in input the iteration count t, creates files for every state of the matrix T and returns 0 if some error occured, 1 otherwise
    */
    int compute(std::size_t t);

};

#endif