#ifndef HH_LATTICE_MATRIX
#define HH_LATTICE_MATRIX

#include <iostream>
#include <string>
#include <vector>

#include "latticeNode.hpp"

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
    const std::string input_file_path;
    
    /**
     * Path to the output directory where all the files will be saved
    */
    const std::string output_dir_path;
    
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
     * Saves output data to the file passed as a constructor to the lattice object
    */
    virtual void save_output_data() const = 0; 
    
    /**
     * Reads the input file producing the lattice.
    */
    virtual void read_input_file() = 0;

    /**
     * Logs specific data regarding the lattice
     * - Number of elements
    */
    virtual void log_specific_data() const = 0;
public:
    
    /**
     * @param input_file_path path to the input file, containing the representation of the simulation domain, used to produce the lattice
     * @param output_file_path_ path to where the output files will be stored
     * @param dimensions_ number of dimensions of the lattice
     * @param velocity_set_ the velocity set used in the lattice
    */
    Lattice(const std::string& input_file_path_, const std::string& output_dir_path_, const int dimensions_, const VelocitySet& velocity_set_);

    virtual ~Lattice() = default;

    /**
     * Performs a step of the simulation of the fluid system by applying the 
     * Lattice Boltzmann discretized equations for each node of the lattice. 
    */
    virtual void perform_simulation_step() = 0;
    /**
     * Performs the initialization of the lattice nodes.
    */
    virtual void initialize_lattice() = 0;

    /**
     * Logs data regarding the lattice: dimensions, path to files and directories,  
    */
    void log_data() const;
};

#endif