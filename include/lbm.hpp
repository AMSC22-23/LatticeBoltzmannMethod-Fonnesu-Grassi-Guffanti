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

// TODO: ADD VISCOSITY V = Cs^2*(tau-dt/2)
// TODO: ADD FUNCTION OBJECT TO MANAGE INLET FIELDS
// TODO: ADD OUTPUT FREQUENCY
// TODO: ADD REYNOLDS NUMBER

class lbm
{
private:
    
    /**
     * The chosen velocity set
    */
    VelocitySet velocity_set;

    /**
     * Lattice of the simulation
    */
    std::unique_ptr<Lattice> lattice_ptr;

    /**
     * Collision model used
    */
    std::shared_ptr<CollisionModel> collision_ptr;

    /**
     * Boundary model used to calculate collisions at the boundary
    */
    std::shared_ptr<Boundary> boundary_ptr;
    
    /**
     * Time constant
    */
    double const tau;

    /**
     * Time lapse
    */
    double const delta_t;

    /*
        function called by constructor setting all inital state
    */
    void load_initial_state();

public:
    /*
        Costructor taking input data from file
    */
    lbm(size_t D, int Q, const std::string& input_dir_path, const std::string& collision_model, const std::string& boundary_model);

    /*
        Solving function, taked in input the iteration count t, creates files for every state of the matrix T and returns 0 if some error occured, 1 otherwise
    */
    int compute(const std::size_t n_iter);

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

    void print()
    {
        lattice_ptr->save_output_data(0);
    }

};

#endif