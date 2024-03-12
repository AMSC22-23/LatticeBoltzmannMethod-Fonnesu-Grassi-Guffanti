#ifndef LLALBM_LATTICEREADER_HPP
#define LLALBM_LATTICEREADER_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
#include <fstream>
#include <algorithm>
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
     * @brief Given a header line as a unique string, this method splits the line into multiple substrings
     * based on a delimiter
     * 
     * @param header Header line
     * @param delimiter Delimiter of the single words
     * @returns Vector of splitted words 
     */
    std::vector<std::string> split_header_line(const std::string& header, const std::string& delimiter)
    {
        std::vector<std::string> splitted;


        std::size_t prev_pos = 0;
        std::size_t pos = header.find(delimiter, prev_pos);

        while (pos != std::string::npos)
        {
            splitted.push_back(header.substr(prev_pos, pos-prev_pos));
            prev_pos = pos + delimiter.size();
            pos = header.find(delimiter, prev_pos);
        }

        splitted.push_back(header.substr(prev_pos, header.size() - 1));
        return splitted;
    }

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
        assert(in.is_open() && "ERROR, the input file could not be opened.");

        // If the file has been opened correctly, then the header is checked.
        // A correct Matrix Market Header is as follows

        logger.info("File was opened, validating header");
        
        std::string header;
        std::stringstream string_stream;
        
        // Check consistency of the header 
        //  - template dimensions == header dimensions
        //  - extensions are positive
        //  - sum of elements is equal to product of extensions
        std::size_t file_dimensions;
        

        // The first element should be the number of dimensions
        std::getline(in, header);
        string_stream << header;
        string_stream >> file_dimensions;
        string_stream.clear();

        logger.info("Physical Dimensions: " + std::to_string(file_dimensions));
        assert(dim == file_dimensions && "ERROR: File provided dimensions and template dimensions do not match."); 

        // Then there should be the list of extensions of the dimensions.
        std::getline(in, header);
        logger.info("Extension of each dimension: " + header);
        // Call the split_header_line method to extract the extensions of the dimensions
        std::vector<std::string> extensions = split_header_line(header, " ");
        assert(extensions.size() == dim && "ERROR: Provided extensions do not match dimensions in number.");

        // Print each element of the extensions vector
        std::size_t total_elems = 1;

        for (std::size_t i = 0; i  < extensions.size(); i++)
        {
            string_stream << extensions[i];
            string_stream >> lattice_dimensions[i];
            
            assert(lattice_dimensions[i] > 0 && "ERROR: Extension of the dimensions must be positive.");

            total_elems = total_elems * lattice_dimensions[i];

            string_stream.clear();
        }
        logger.info("Lattice dimensions stored.");

        // Now read the next line of the header file, which must contain the number of elements of each node type.
        // As previously stated, the number of different node types is 5
        // - Fluid nodes
        // - Solid nodes
        // - Boundary nodes
        // - Inlet nodes
        // - Outlet nodes
        // - Obstacle nodes
        // => refer to the InputNodeTypes enum for the order of the types.
        
        // Load the data in a std::vector and split it
        std::getline(in, header);
        std::vector<std::string> node_types = split_header_line(header, " ");
        std::vector<std::size_t> node_counts(node_types.size());
        
        std::size_t all_different_nodes = 0;
        logger.info("Nodes of different types: " + header);
        for (std::size_t i = 0; i < node_types.size(); ++i)
        {
            string_stream << node_types[i];
            string_stream >> node_counts[i];
            string_stream.clear();
            
            assert(node_counts[i] >= 0 && "Error: Number of node types must be non-negative.");
            all_different_nodes += node_counts[i];
        }
        assert(all_different_nodes == total_elems && "Error: Total number of elements and described number of nodes do not match");
        logger.info("==Header read correcly==");

        // Data structures allocation
        // Allocate space for intermediate tensor of types

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