#ifndef HH_LATTICE_2D
#define HH_LATTICE_2D

#include <cassert>
#include <fstream>
#include <vector>
#include <filesystem>
#include <numeric>
#include <cmath>
#include <string>
#include <algorithm>
#include <execution>
#include <omp.h>

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

    static constexpr double one_third = 1.0/3.0;
    
    std::size_t num_of_threads = 1;
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

    virtual void calculate_equilibrium() override;

    /**
     * Executes streaming of the populations
    */
    virtual void perform_streaming() override;

    virtual void perform_boundary_collisions() override;

    virtual void calculate_macroscopic_quantities() override;
   
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
     * @param boundary_model model used to calculate collisions on boundaries
     * @param reynolds_ reynolds constant
    */ 
    Lattice2D(const std::string& input_file_path_,
        const std::string& output_dir_path_,
        const VelocitySet& velocity_set_,
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

    /**
     * Sets the number of threads
    */
    std::size_t& set_omp_num_threads();

    /**
     * returns the dimension of the lattice as [width, height]
    */
    virtual const std::vector<std::size_t> get_lattice_dimension() const override;

    /**
      * Attaches a collision model to be used in the lattice
    */
    virtual void attach_collision_model(const std::shared_ptr<CollisionModel>& model) override;
};

#endif // HH_LATTICE_2D