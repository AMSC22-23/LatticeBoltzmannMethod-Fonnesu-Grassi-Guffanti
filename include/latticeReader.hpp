#ifndef HH_LATTICE_READER
#define HH_LATTICE_READER

#include <iostream>
#include <string>
#include <filesystem>
#include <cassert>
#include <string>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

#include "latticeNode.hpp"
#include "aliases.hpp"
#include "boundary.hpp"

/**
 * Abstract class that virtual methods that are overridden to implement the 
 * reading of a lattice structure and the values used to initialize the lattice
 * @author Luca Guffanti
*/
class LatticeReader
{
protected:
    /**
     * Path to the input directory which contains all the input files
    */
    const std::string input_dir_path;
    /**
     * Path to the input matrix describing the lattice structure
    */
    const std::string input_lattice_path;
    /**
    *  Path to the input matrix describing the lattice density field
    */
    const std::string input_rho_path;
    /**
    *  Path to the input matrix describing the lattice velocity fields 
    */
    const std::string input_u_path;
    
private:
    /**
     * Method that validates the path passed as input
    */
    bool validate_path();
public:
    LatticeReader(const std::string& input_dir_path_);
    ~LatticeReader() = default;
};

/**
 * Implementation of the abstract LatticeReader class
 * used to read the sets of matrices that describe a 2D lattice and
 * the initialization values
*/
class LatticeReader2D : public LatticeReader
{
private:
    bool load_2D_velocity_matrix(LatticeGrid2D& lattice, const std::string& path, const std::size_t& dir);
public:
    /**
     * Constructs the object
     * @param input_dir_path the path to the input directory
    */
    LatticeReader2D(const std::string& input_dir_path_);
    ~LatticeReader2D() = default;

    /**
     * Reads the lattice node structure and populates
     * @param lattice the lattice to be constructed
     * @param boundary_list the list of boundary nodes that will be populated and used to manage boundary collisions.
     * @param width of the lattice
     * @param height of the lattice
    */
    bool read_lattice_structure(LatticeGrid2D& lattice, BoundaryList2D& boundary_list,std::size_t& width, std::size_t& height);
    /**
     * Reads the input density field from file and coherently updates lattice nodes.
     * @param lattice lattice containing the nodes whose density is to be updated
    */
    bool read_lattice_input_rho(LatticeGrid2D& lattice);
    /**
     * Reads the input velocity field from file and coherently updates lattice nodes.
     * @param lattice lattice containing the nodes whose velocities are to be updated
    */
    bool read_lattice_input_velocities(LatticeGrid2D& lattice);
};


#endif // HH_LATTICE_READER