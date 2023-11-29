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
#include "boundary.hpp"

class lbm
{
private:
    
    VelocitySet velocity_set;

    std::unique_ptr<Lattice> lattice_ptr;
    std::unique_ptr<CollisionModel> collision_ptr;
    std::unique_ptr<Boundary> boundary_ptr;
    // TODO: aggiungere tutte le costanti del caso
    double const tau;
    double const delta_t;



    /*
        function called by constructor setting all inital state
    */
    void load_initial_state();

public:
    /*
        Costructor taking input data from file
    */
    lbm(size_t D, int Q, const std::string& input_file_path, const std::string& collision_model, const std::string& boundary_model);

    /*
        Solving function, taked in input the iteration count t, creates files for every state of the matrix T and returns 0 if some error occured, 1 otherwise
    */
    int compute(const int n_iter);

    double const get_tau() const{
        return tau;
    }

    double const get_delta_t() const{
        return delta_t;
    }

    VelocitySet const get_set() const{
        return velocity_set;
    }

    /*std::unique_ptr<Lattice>& set_collision_ptr(){
        return lattice_ptr;
    }

    std::unique_ptr<CollisionModel>& set_collision_ptr(){
        return collision_ptr;
    }

    std::unique_ptr<Boundary>& set_collision_ptr(){
        return boundary_ptr;
    }*/

};

#endif