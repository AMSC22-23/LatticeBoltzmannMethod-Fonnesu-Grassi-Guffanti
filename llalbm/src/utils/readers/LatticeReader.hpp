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
    std::vector<std::string> split_file_line(const std::string& line, const std::string& delimiter)
    {
        std::vector<std::string> splitted;


        std::size_t prev_pos = 0;
        std::size_t pos = line.find(delimiter, prev_pos);

        while (pos != std::string::npos)
        {
            splitted.push_back(line.substr(prev_pos, pos-prev_pos));
            prev_pos = pos + delimiter.size();
            pos = line.find(delimiter, prev_pos);
        }

        splitted.push_back(line.substr(prev_pos, line.size() - 1));
        return splitted;
    }


// =======================================================================================================
//                                  TEMPLATIED STRUCTS FOR TENSORS
// =======================================================================================================

    /**
     * @brief Class that iterates through a tensor and populates vectors of coordinates
     * based on the value of the tensor at the given coordinates, which described a type of node.
     * @note The class is templated on the type of the tensor, the number of ranks and the current dimension
     * in order to allow for a generic implementation of the tensor iteration, especially independent of the
     * number of dimensions of the tensor.
     * 
     * This choice allows for a more flexible implementation of the iteration, which is substistuted by a
     * recursive call to the same method, with the current dimension incremented by one.
     * Last but not least, the class is specialized for the case in which the current dimension is equal to the
     * number of ranks, in which case the method populates the vectors with the coordinates of the nodes, passed
     * as a parameter and built during the sequence of recursive calls.
     * 
     * @tparam TensorType 
     * @tparam ranks 
     * @tparam current_dim 
     */
    template<typename TensorType, std::size_t ranks, std::size_t current_dim = 0>
    struct TensorLoopAccessor
    {
        /**
         * @brief Method that iterates through a specific dimension
         * of the tensor, described by the template parameter, calling itself recursively to 
         * isolate, at the end, a single element.
         * 
         * @param tensor Tensor to be iterated
         * @param fluid_nodes_coord Buffer containing the coordinates of fluid nodes
         * @param boundary_coord Buffer containing the coordinates and types of boundary nodes
         * @param inlet_nodes_coord Buffer containing the coordinates and types of inlet nodes
         * @param outlet_nodes_coord Buffer containing the coordinates and types of outlet nodes
         * @param obstacle_nodes_coord Buffer contatining the coordinates and types of obstacle nodes
         * @param indices Array containing the current indices of the tensor
         * @param occupations Array containing the actual dimension of the vectors of each type of node
         * 
         */
        static void vector_population_loop(
            const TensorType& tensor,
            std::vector<Point<ranks>>& fluid_nodes_coord,
            std::vector<BoundaryPoint<ranks>>& boundary_coord,
            std::vector<BoundaryPoint<ranks>>& inlet_nodes_coord,
            std::vector<BoundaryPoint<ranks>>& outlet_nodes_coord,
            std::vector<BoundaryPoint<ranks>>& obstacle_nodes_coord,
            Eigen::Index indices[],
            std::array<std::size_t, 5>& occupations
        )
        {
            auto size = tensor.dimensions()[current_dim];
            for (auto idx = 0; idx < size; ++idx)
            {
                indices[current_dim] = idx;
                TensorLoopAccessor<TensorType, ranks, current_dim + 1>::vector_population_loop(
                    tensor,
                    fluid_nodes_coord,
                    boundary_coord,
                    inlet_nodes_coord,
                    outlet_nodes_coord,
                    obstacle_nodes_coord,
                    indices,
                    occupations
                );
            }
        }
    };

    /**
     * @brief Specialization of the TensorLoopAccessor class for the case in which the current dimension
     * is equal to the number of ranks of the tensor. In this case, the method populates the vectors with the
     * coordinates of the nodes, passed as a parameter and built during the sequence of recursive calls.
     * 
     * @tparam TensorType Type of the tensor
     * @tparam ranks Total number of ranks of the tensor
     */
    template<typename TensorType, std::size_t ranks>
    struct TensorLoopAccessor<TensorType, ranks, ranks>
    {
        /**
         * @brief Method that populates the vectors with the coordinates of the nodes, passed as a parameter and built during the sequence of recursive calls.
         * 
         * @note The method is specialized for the case in which the current dimension is equal to the number of ranks, 
         * which indicates that the tensor has been reduced to a single value.
         * 
         * @param tensor Tensor to be iterated
         * @param fluid_nodes_coord Buffer containing the coordinates of fluid nodes
         * @param boundary_coord Buffer containing the coordinates and types of boundary nodes
         * @param inlet_nodes_coord Buffer containing the coordinates and types of inlet nodes
         * @param outlet_nodes_coord Buffer containing the coordinates and types of outlet nodes
         * @param obstacle_nodes_coord Buffer contatining the coordinates and types of obstacle nodes
         * @param indices Array containing the current indices of the tensor
         * @param occupations Array containing the actual dimension of the vectors of each type of node
         */
        static void vector_population_loop(
            const TensorType& tensor,
            std::vector<Point<ranks>>& fluid_nodes_coord,
            std::vector<BoundaryPoint<ranks>>& boundary_coord,
            std::vector<BoundaryPoint<ranks>>& inlet_nodes_coord,
            std::vector<BoundaryPoint<ranks>>& outlet_nodes_coord,
            std::vector<BoundaryPoint<ranks>>& obstacle_nodes_coord,
            Eigen::Index indices[],
            std::array<std::size_t, 5>& occupations
        )
        {
            // As we reduce the tensor to a single value, we just need to access the element
            // and update the corresponding array.

            std::size_t type = tensor(indices);
            Point<ranks> p;
            BoundaryPoint<ranks> bp;
            for (std::size_t i = 0; i < ranks; ++i)
            {
                p.coords[i] = indices[i];
                bp.coords[i] = indices[i];
            }

            bp.type = btypes::NONE;
            switch(type)
            {
                case InputNodeType::FLUID:
                    fluid_nodes_coord.at(occupations[0]) = p;
                    occupations[0]++;
                    break;
                case InputNodeType::BOUNDARY:
                    boundary_coord.at(occupations[1]) = bp;
                    occupations[1]++;
                    break;
                case InputNodeType::INLET:
                    inlet_nodes_coord.at(occupations[2]) = bp;
                    occupations[2]++;
                    break;
                case InputNodeType::OUTLET:
                    outlet_nodes_coord.at(occupations[3]) = bp;
                    occupations[3]++;
                    break;
                case InputNodeType::OBSTACLE:
                    obstacle_nodes_coord.at(occupations[4]) = bp;
                    occupations[4]++;
                    break;
            }
        }
    };


    /**
     * @brief Wrapper around a single static method that generically loads a tensor from an
     * input file stream. While programmed for the specific case of a simil-mtx file format,
     * this function behaves correctly whenever tensor elements are described as a list of coordinates and the value
     * @note Templating is in the type of tensor, so that the function can be reused
     * 
     * @tparam TensorType Type of the tensor to be loaded
     */
    template<typename TensorType>
    struct TensorLoader
    {
        /**
         * @brief Loads a tensor from an input file stream. Each line of the file (after the header)
         * should contain the appropriate number coordinates and the value of the tensor correspondnting to those
         * coordinates
         * 
         * @param in Input file stream
         * @param tensor Tensor to be loaded
         * @param total_elems Total number of expected lines
         * @param dims Number of ranks of the tenso
         * @param coords Vector, of type @code
         */
        static void load_tensor_from_file(std::ifstream& in, TensorType& tensor, const std::size_t& total_elems, const std::size_t& dims, Eigen::Index coords[])
        {
            std::string file_line;
            std::vector<std::string> splitted;
            std::stringstream string_stream;

            std::size_t type;

            std::size_t count = 0;
            assert(in.is_open() && "ERROR, the input file could not be opened.");
            
            while(count < total_elems)
            {
                getline(in, file_line);
                splitted = split_file_line(file_line, " ");
                
                // Load the dimensions in the array
                for (std::size_t i = 0; i < dims; ++i)
                {
                    string_stream << splitted[i];
                    string_stream >> coords[i];
                    string_stream.clear();
                }
                
                string_stream << splitted[dims];
                string_stream >> type;
                string_stream.clear();

                tensor(coords) = type;
                ++count;
            }
        }
    };
// =======================================================================================================
// =======================================================================================================
// =======================================================================================================

// =======================================================================================================
//                          POST-LOADING NODE IDENTIFICATION
// =======================================================================================================

    /**
     * @brief Checks wheter the vector of nodes passed as a parameter contains nodes that are at the bounds of the
     * lattice. 
     * @note Whilts available to all types of nodes, this method is thought to be used with inlet and outlet nodes:
     * a lattice that presents inlets and outlets which are not at the physical extremes of the domain is considered 
     * incorrect 
     * 
     * @tparam dim Physical dimensions of the lattice 
     * @param lattice_dimensions Extension of each lattice dimension
     * @param nodes Vector of nodes to be checked
     */
    template<std::size_t dim>
    void verify_lattice_bounds(const std::array<Eigen::Index, dim>& lattice_dimensions, const std::vector<BoundaryPoint<dim>>& nodes)
    {
        // For each node, check that each of the coordinate is either at the beginning or at the end of the corresponding lattice dimension.
        bool verified;
        const std::size_t size = nodes.size();
        for (std::size_t i = 0; i < size; ++i)
        {
            verified =false;
            for (std::size_t j = 0; j < dim && !verified; ++j)
            {
                verified = nodes[i].coords[j] == 0 || nodes[i].coords[j] == lattice_dimensions[j] - 1;
            }
            assert(verified && "ERROR: A node is not at the boundary while it must be");
        }
    }

    /**
     * @brief Identifies the boundary type of each node based on its coordinates in a 2D lattice
     * @note The method is used after the vectors have been fully populated, with the boundary type of 
     * node being the only thing that is to be identified.
     * 
     * @param lattice_dimensions Extension of each lattice dimension
     * @param nodes Vector of nodes to be checked
     */
    void identify_node_type(const std::array<Eigen::Index, 2>& lattice_dimensions, std::vector<BoundaryPoint<2>>& nodes)
    {
        // For each node, check the position of the node and assign the type accordingly. 
        // In the identification, it's to be noted that the lattice dimensions proceed as follows (just like when describing a matrix):
        //  - The first dimension is the row index
        //  - The second dimension is the column index

        std::size_t size = nodes.size();
        for (std::size_t i = 0; i < size; ++i)
        {
            // The node is at the upper boundary if the row index is 0
            if (nodes[i].coords[0] == 0)
            {
                // Then, if the column index is 0, the node is at the top left corner
                if (nodes[i].coords[1] == 0)
                {
                    nodes[i].type = btypes::LEFT_TOP_C;
                }
                // If the column index is equal to the dimension - 1, the node is at the top right corner
                else if (nodes[i].coords[1] == lattice_dimensions[1] - 1)
                {
                    nodes[i].type = btypes::RIGHT_TOP_C;
                }
                // Otherwise, the node is at the top wall
                else
                {
                    nodes[i].type = btypes::TOP_WALL;
                }
            }
            // The node is at the lower boundary if the row index is equal to the dimension - 1
            else if (nodes[i].coords[0] == lattice_dimensions[0] - 1)
            {
                // Then, if the column index is 0, the node is at the bottom left corner
                if (nodes[i].coords[1] == 0)
                {
                    nodes[i].type = btypes::LEFT_BOTTOM_C;
                }
                // If the column index is equal to the dimension - 1, the node is at the bottom right corner
                else if (nodes[i].coords[1] == lattice_dimensions[1] - 1)
                {
                    nodes[i].type = btypes::RIGHT_BOTTOM_C;
                }
                // Otherwise, the node is at the bottom wall
                else
                {
                    nodes[i].type = btypes::BOTTOM_WALL;
                }
            }
            // Now, the node being a wall or a corner has already been checked, so we must just check whether the
            // node is on the left or on the right of the lattice.
            // The node is at the left wall if the column index is 0 and the row index is not 0 or the dimension - 1
            else if (nodes[i].coords[1] == 0 && nodes[i].coords[0] != 0 && nodes[i].coords[0] != lattice_dimensions[0] - 1)
            {
                nodes[i].type = btypes::LEFT_WALL;
            }
            // Finally, the node is at the right wall if the column index is equal to the dimension - 1 and the row index is not 0 or the dimension - 1
            else if (nodes[i].coords[1] == lattice_dimensions[1] - 1 && nodes[i].coords[0] != 0 && nodes[i].coords[0] != lattice_dimensions[0] - 1)
            {
                nodes[i].type = btypes::RIGHT_WALL;
            }
        }
        
            
    }


// =======================================================================================================
//                          MAIN NAMESPACE FUNCTION
// =======================================================================================================

    /**
     * @brief Method that reads lattice data from an input file that is 
     * stored in matrix market format, describing the lattice nodes as in a matrix.
     * As the systems allows the construction of tensors, the code must behave in
     * a generic way. Vectors containing coordinates of data in the lattice will be 
     * instantiated in the method. 
     * 
     * @note To operate correctly the file should be structured in the following way.
     * 
     * NUMBER OF PHYSICAL DIMENSIONS
     * SPACE SEPARATED EXTENSIONS OF EACH DIMENSION, PROCEEDING AS FOLLOWS (JUST LIKE WHEN DEFINING A MATRIX):
     * - first number  : height of the lattice
     * - second number : width of the lattice
     * - third number  : depth of the lattice
     * ... and so on
     * 
     * NUMBER OF LATTICE NODES PER TYPE
     * 
     * SPACE SEPARATED COORDINATES OF THE NODES AND TYPE, PROCEEDING AS FOLLOWS (JUST LIKE WHEN POPULATING A MATRIX):
     * - first number  : row index 
     * - second number : column index
     * - third number  : depth index
     * ... and so on
     * 
     * For instance, considering only three types of nodes, a file
     * 3
     * 4 3 5
     * 1 10 30
     * ...
     * 
     * is not valid has the total number of points does not reach the product
     * of the extensions of the dimensions.
     * 
     * While a file
     * 2
     * 10 20
     * 50 50 100
     * 1 2 1
     * 3 6 3
     * 7 19 4
     * ...
     * 
     * is valid as the total number of points is equal to the product of the extensions, and the coordinates indicate
     * - a solid node at row 1 and column 2
     * - an inlet node at row 3 and column 6
     * - an outlet node at row 7 and column 19
     * 
     * Nonetheless, values referring to FLUID NODES are not indicated: the leading
     * assumption is that the majority of nodes will be fluid, so to compress the input
     * file the deduction of where a fluid node is located is done at runtime.
     * This choice reduces the spatial occupancy of the file, which should not be cumbersome,
     * but in turn increases the time needed to build the data structures.
     * 
     * @tparam dim Number of physical dimensions
     * @param path Path to the input file
     * @param fluid_nodes_coord Buffer containing the coordinates of fluid nodes
     * @param boundary_coord Buffer containing the coordinates and types of boundary nodes
     * @param inlet_nodes_coord Buffer containing the coordinates and types of inlet nodes
     * @param outlet_nodes_coord Buffer containing the coordinates and types of outlet nodes
     * @param obstacle_nodes_coord Buffer contatining the coordinates and types of obstacle nodes
     * @param lattice_dimensions Number of element per dimension
     */
    template<std::size_t dim>
    void read_lattice_file(
        const std::string& path,
        std::vector<Point<dim>>& fluid_nodes_coord,
        std::vector<BoundaryPoint<dim>>& boundary_coord,
        std::vector<BoundaryPoint<dim>>& inlet_nodes_coord,
        std::vector<BoundaryPoint<dim>>& outlet_nodes_coord,
        std::vector<BoundaryPoint<dim>>& obstacle_nodes_coord,
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
        // Call the split_file_line method to extract the extensions of the dimensions
        std::vector<std::string> extensions = split_file_line(header, " ");
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
        std::vector<std::string> node_types = split_file_line(header, " ");
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
        Eigen::Tensor<std::size_t, dim> intermediate_types(lattice_dimensions);

        // Allocate space for the single vectors
        fluid_nodes_coord.resize(node_counts[InputNodeType::FLUID]);
        boundary_coord.resize(node_counts[InputNodeType::BOUNDARY]);
        inlet_nodes_coord.resize(node_counts[InputNodeType::INLET]);
        outlet_nodes_coord.resize(node_counts[InputNodeType::OUTLET]);
        obstacle_nodes_coord.resize(node_counts[InputNodeType::OBSTACLE]);

        logger.info("Fluid nodes vector resized to " + std::to_string(node_counts[InputNodeType::FLUID]));
        logger.info("Boundary nodes vector resized to " + std::to_string(node_counts[InputNodeType::BOUNDARY]));
        logger.info("Inlet nodes vector resized to " + std::to_string(node_counts[InputNodeType::INLET]));
        logger.info("Outlet nodes vector resized to " + std::to_string(node_counts[InputNodeType::OUTLET]));
        logger.info("Obstacles nodes vector resized to " + std::to_string(node_counts[InputNodeType::OBSTACLE]));

        // Load the contents of the file in memory
        std::size_t non_fluids = all_different_nodes - fluid_nodes_coord.size();
        
        Eigen::Index indices[dim];
        std::array<std::size_t, 5> occupations;
        std::fill(occupations.begin(), occupations.end(), 0);
        TensorLoader<Eigen::Tensor<std::size_t, dim>>::load_tensor_from_file(in, intermediate_types, non_fluids, file_dimensions, indices);
        // Then, populate with the correct coordinates the vectors.

        logger.info("==Read completed==");
        in.close();

        // Fluid node identification
        // Iterating through the intermediate tensor, identify the fluid
        // nodes based on their id, and consenquently add the coordinates to the 
        // fluid node vector


        TensorLoopAccessor<Eigen::Tensor<std::size_t, dim>, dim>::vector_population_loop(
            intermediate_types,
            fluid_nodes_coord,
            boundary_coord,
            inlet_nodes_coord,
            outlet_nodes_coord,
            obstacle_nodes_coord,
            indices,
            occupations
        );


        // Element analysis
        // Analyze only the inlet nodes:
        //  - Check that they are at the extremes of the lattice, otherwise the lattice cannot be
        //      considered usable 
        //  - For the inlets and outlets, only consider top, bottom, right and left walls

        // Inlet nodes
        verify_lattice_bounds<dim>(lattice_dimensions, inlet_nodes_coord);
        // Outlet nodes
        verify_lattice_bounds<dim>(lattice_dimensions, outlet_nodes_coord);

        // Identify the type of inlet and outlet nodes
        identify_node_type(lattice_dimensions, inlet_nodes_coord);
        identify_node_type(lattice_dimensions, outlet_nodes_coord);
        logger.info("==Nodes identified properly==");
        logger.info("LATTICE READING COMPLETED");
    }

// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
    

};

#endif // LLALBM_LATTICEREADER_HPP