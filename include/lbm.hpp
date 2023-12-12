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
    double tau;

    /**
     * Time interval
    */
    double delta_t;

    /**
     * Frequency at which the output will be presented
    */
    const std::size_t frequency;

    /**
     * Reynolds number
    */
    const double re;

    /**
    * Function called by constructor setting all inital state
    */
    void load_initial_state();

public:
    /**
    * Costructor taking input data from file
    * @param D the number of dimensions
    * @param reynolds_ the reynolds number
    * @param input_dir_path directory to input containing matrices
    * @param collision_model collision model
    * @param frequency_ frequency at which the output will be presented
    */
    lbm(const std::size_t& D,
        const double reynolds_, 
        const std::string& input_dir_path, 
        const std::string& collision_model,
        const std::size_t& frequency_);
    /*
        Solving function, taked in input the iteration count t, creates files for every state of the matrix T and returns 0 if some error occured, 1 otherwise
    */
    void compute(const double time);

    double const get_tau() const
    {
        return tau;
    }

    double const get_delta_t() const
    {
        return delta_t;
    }

    VelocitySet const get_set() const
    {
        return velocity_set;
    }
};

#endif