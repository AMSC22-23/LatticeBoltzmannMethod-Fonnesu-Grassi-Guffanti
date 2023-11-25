#ifndef HH_LATTICE_2D
#define HH_LATTICE_2D

#include <cassert>
#include <fstream>
#include <vector>
#include <string>

#include "lattice.hpp"
#include "latticeNode.hpp"

class Lattice2D : public Lattice
{
private:
    /**
     * Dimensions of the lattice
    */
    static constexpr std::size_t dim = 2;
   
    /**
     * Lattice, represented as a two-dimensional array
    */
    std::vector<std::vector<LatticeNode<dim>>> lattice;
   
    /**
     * Saves the output data to file dor the 2D lattice
     * @note Implementation of Lattice class equivalent virtual method
    */
    virtual void save_output_data() const override; 
   
    /**
     * Reads an input file describing a 2D domain
     * @note Implementation of Lattice class equivalent virtual method
    */
    virtual void read_input_file() override;
   
    /**
     * Logs specific data regarding the lattice
     * - Number of elements per dimension
    */
    virtual void log_specific_data() const override;
public:
    
    /**
     * @param input_file_path path to the input file
     * @param output_dir_path path to the output directory
    */
    Lattice2D(const std::string& input_file_path_, const std::string& output_dir_path_);
    virtual ~Lattice2D() = default;

    /**
     * Performs a step of the simulation of the fluid system by applying the 
     * Lattice Boltzmann discretized equations for each node of the 2D lattice. 
    */
    virtual void perform_simulation_step() override;
    
    /**
     * Initializes the 2D lattice to default values of populations and fields.
    */
    virtual void initialize_lattice() override;
};

#endif // HH_LATTICE_2D