#ifndef HH_LATTICE_MATRIX
#define HH_LATTICE_MATRIX

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <cassert>

#include "latticeNode.hpp"
#include "boundary.hpp"
#include "inletInitializer.hpp"

/**
 * Abstract representation of the lattice of points, realized into two
 * different implementations for the 2D and the 3D use-case.
 * @author Luca Guffanti
*/
class Lattice
{
protected:

    /**
     * Path to the input file where the domain of the computation is described
    */
    const std::string input_dir_path;
    
    /**
     * Path to the output directory where all the files will be saved
    */
    std::string output_dir_path;
    
    /**
     * Number of dimensions of the lattice
    */
    const std::size_t dimensions;
    
    /**
     * Name of the domain
    */
    std::string lattice_name;
    
    /**
     * The chosen velocity set
    */
    const VelocitySet velocity_set;

    /**
     * Collision model used
    */
    std::shared_ptr<CollisionModel> collision_model;

    /**
     * Time constant
    */
    double tau;

    /**
     * Time lapse
    */
    double delta_t;

    /**
     * Logs specific data regarding the lattice
     * - Number of elements
    */
    virtual void log_specific_data() const = 0;

    /**
     * Executes streaming of the populations
    */
    virtual void perform_streaming() = 0;

    /**
    * Creates the directory in which the results will be stored
    */
    void create_output_directory();
public:
    
    /**
     * @param input_dir_path path to the input file, containing the representation of the simulation domain, used to produce the lattice
     * @param output_file_path_ path to where the output files will be stored
     * @param dimensions_ number of dimensions of the lattice
     * @param velocity_set_ the velocity set used in the lattice
     * @param collision_model model used to calculate collisions
     * @param reynolds_ reynolds constant
    */
    Lattice(const std::string& input_dir_path_, 
        const std::string& output_dir_path_, 
        const int dimensions_, 
        const VelocitySet& velocity_set_, 
        std::shared_ptr<CollisionModel> collision_model_, 
        const double reynolds_);

    virtual ~Lattice() = default;

    /**
     * Performs a step of the simulation of the fluid system by applying the 
     * Lattice Boltzmann discretized equations for each node of the lattice. 
    */
    virtual void perform_simulation_step() = 0;
    
    /**
     * Saves output data to the file passed as a constructor to the lattice object
     * @param iteration_count the current iteration number, in order to produce the correct output file
    */
    virtual void save_output_data(std::size_t iteration_count) const = 0; 

    /**
     * Performs the initialization of the lattice nodes.
    */
    virtual void initialize_lattice() = 0;

    /**
     * Sets the inlet fields of the lattice;
    */
    virtual void set_inlets(const std::size_t iterations) = 0;

    /**
     * Logs data regarding the lattice: dimensions, path to files and directories,  
    */
    void log_data() const;
};

#endif