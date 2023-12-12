#ifndef HH_LATTICE_2D
#define HH_LATTICE_2D

#include <cassert>
#include <fstream>
#include <vector>
#include <filesystem>
#include <numeric>
#include <cmath>
#include <string>

#include "aliases.hpp"
#include "lattice.hpp"
#include "latticeNode.hpp"
#include "collisionModel.hpp"
#include "boundary.hpp"
#include "latticeReader.hpp"
#include "inletInitializer.hpp"

/**
 * 2D implementation of the Lattice abstract class
 * @author Luca Guffanti
*/
class Lattice2D : public Lattice
{
private:

    /**
     * Dimensions of the lattice
    */
    static constexpr std::size_t dim = 2;
    
    /**
     * Boundary model
    */
    NEBB boundary_model;
    /**
     * Width of the lattice
    */
    std::size_t lattice_width;
    /**
     * Height of the lattice
    */
    std::size_t lattice_height;
   
    /**
     * Lattice, represented as a two-dimensional array
    */
    LatticeGrid2D lattice;

    /**
     * Object that implements methods to read a 2D lattice and various initialization values
    */
    std::unique_ptr<LatticeReader2D> lattice_reader;

    /**
     * Object that initializes lattice inlet fields
    */
    std::unique_ptr<InletInitializer<2>> inlet_initializer;

    /**
     * Vector of tuples of boundaries and open boundaries and their coordinates
    */
    BoundaryList2D boundary_list;

    /**
     * Executes streaming of the populations
    */
    virtual void perform_streaming() override;
   
    /**
     * Logs specific data regarding the lattice
     * - Number of elements per dimension
    */
    virtual void log_specific_data() const override;

public:
    
    /**
     * @param input_file_path path to the input file
     * @param output_dir_path path to the output directory
     * @param velocity_set_ velocity set used to perform operations on the lattice
     * @param collision_model model used to calculate collisions
     * @param boundary_model model used to calculate collisions on boundaries
     * @param reynolds_ reynolds constant
    */ 
    Lattice2D(const std::string& input_file_path_,
        const std::string& output_dir_path_,
        const VelocitySet& velocity_set_,
        std::shared_ptr<CollisionModel> collision_model_, 
        const double reynolds_);
    virtual ~Lattice2D() = default;

    /**
     * Performs a step of the simulation of the fluid system by applying the 
     * Lattice Boltzmann discretized equations for each node of the 2D lattice. 
    */
    virtual void perform_simulation_step() override;
    
    /**
     * Initializes the 2D lattice to default values of populations and macroscopic variables.
    */
    virtual void initialize_lattice() override;

    /**
     * Saves the output data to file dor the 2D lattice
     * @note Implementation of Lattice class equivalent virtual method
    */
    virtual void save_output_data(std::size_t iteration_count) const override; 

    /**
     * Sets the inlets fields
    */
    virtual void set_inlets(const std::size_t iterations) override;
};

#endif // HH_LATTICE_2D