#ifndef LLALBM_LATTICEREADER_HPP
#define LLALBM_LATTICEREADER_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
#include <fstream>
// ======================================

// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "../../utils/loggers/Logger.hpp"
#include "../../utils/aliases.hpp" 
// =======================================


// // ========= COORDINATES OF ALL NODES =========

// /// @brief List of coordinates of fluid nodes
// std::vector<Point<dim>> fluid_nodes;

// /// @brief List of coordinates of boundary nodes
// std::vector<boundaryPoint<dim>> boundary_coord;

// /// @brief List of coordinates of inlet nodes.
// std::vector<boundaryPoint<dim>> inlet_nodes_coord;

// /// @brief List of coordinates of outlet nodes.
// std::vector<boundaryPoint<dim>> outlet_nodes_coord;

// /// @brief List of coordinates of the edges of obstacle nodes.
// std::vector<boundaryPoint<dim>> obstacle_nodes;

namespace llalbm::util::reader 
{
    using namespace llalbm::util::logger;

    static Logger logger("LatticeReader", std::cout);
    /**
     * @brief Method that reads lattice data from an input file that is 
     * stored in matrix market format, describing the lattice nodes as in a matrix.
     * As the systems allows the construction of tensors, the code must behave in
     * a generic way.\n 
     * 
     * @note To operate correctly the file should be structured in the following way.
     * 
     * NUMBER OF PHYSICAL DIMENSIONS
     * SPACE SEPARATED EXTENSIONS OF EACH DIMENSION
     * NUMBER OF LATTICE NODES PER TYPE
     * SPACE SPEARATEDCOORDINATES OF THE NODES AND TYPE
     * 
     * 
     * For instance a file
     * 3
     * 4 3 5
     * 1 10 30
     * ...
     * 
     * is not valid has the total number of points does not reach the product
     * of the extensions of the dimensions.
     * Nonetheless, values referring to FLUID NODES are not indicated: the leading
     * assumption is that the majority of nodes will be fluid, so to compress the input
     * file the deduction of where a fluid node is located is done at runtime.
     * This choice reduces the spatial occupancy of the file, which should not be cumbersome,
     * but in turn increases the time needed to build the data structures.
     * 
     * @tparam dim Number of physical dimensions
     * @param path Path to the input file
     * @param fluid_nodes Buffer containing the coordinates of fluid nodes
     * @param boundary_coord Buffer containing the coordinates and types of boundary nodes
     * @param inlet_nodes_coord Buffer containing the coordinates and types of inlet nodes
     * @param outlet_nodes_coord Buffer containing the coordinates and types of outlet nodes
     * @param obstacle_nodes Buffer contatining the coordinates and types of obstacle nodes
     * @param lattice_dimensions Number of element per dimension
     */
    template<std::size_t dim>
    void read_lattice_file(
        const std::string& path,
        std::vector<Point<dim>>& fluid_nodes,
        std::vector<boundaryPoint<dim>>& boundary_coord,
        std::vector<boundaryPoint<dim>>& inlet_nodes_coord,
        std::vector<boundaryPoint<dim>>& outlet_nodes_coord,
        std::vector<boundaryPoint<dim>>& obstacle_nodes,
        std::array<Eigen::Index, dim>& lattice_dimensions)
    {
        // In order to read the file, instantiate an input file stream.
        std::ifstream in(path);
        assert(in.is_open() && "ERROR, the LTC input file could not be opened.");

        // If the file has been opened correctly, then the header is checked.
        // A correct Matrix Market Header is as follows

        logger.info(".ltc file was opened, validating header");
        
        std::string header;
        std::stringstream string_stream;
        
        // Check consistency of the header 
        //  - template dimensions == header dimensions
        //  - extensions are positive
        //  - sum of elements is equal to product of extensions
        std::getline(in, header);
        
        // The first element should be the number of dimensions
        logger.info(header);

        // Data structures allocation
        // TODO: Allocate space for intermediate tensor of types
        // TODO: Allocate space for the vectors
        // TODO: Store the extension of each dimension in the lattice_dimensions array

        // File reading
        // TODO: Read the file line by line, extract coordinates and file type.
        // Based on the file type add an element to the correct vector.
        
        // Fluid node identification
        // TODO: Iterating through the intermediate tensor, identify the fluid
        // nodes based on their id, and consenquently add the coordinates to the 
        // fluid node vector

        // Element analysis
        // TODO: Analyze each element of each non-fluid vector and identify the type
        // of element.

    }


};

#endif // LLALBM_LATTICEREADER_HPP