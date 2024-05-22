#ifndef LLALBM_NODEIDENTIFICATION_HPP
#define LLALBM_NODEIDENTIFICATION_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
#include <fstream>
// ====================================
// =========== EIGEN INCLUDES ==========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Core>
#include <Eigen/Dense>
// ====================================
// =========== LLALBM INCLUDES ==========
#include "../aliases.hpp"
// ====================================


namespace llalbm::util
{
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

}


#endif