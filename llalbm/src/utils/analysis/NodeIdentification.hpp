#ifndef LLALBM_NODEIDENTIFICATION_HPP
#define LLALBM_NODEIDENTIFICATION_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
#include <functional>
#include <fstream>
#include <map>
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
     * @brief Identifies the boundary type of each node in the most general d-dimensional computational domain.
     * 
     * @tparam d 
     * @param lattice_dimensions 
     * @param nodes 
     */
    template<std::size_t d>
    void identify_node_type(const std::array<Eigen::Index, d>& lattice_dimensions, std::vector<BoundaryPoint<d>>& nodes)
    {
        return;
    }

    /**
     * @brief Identifies the boundary type of each node based on its coordinates in a 2D lattice
     * @note The method is used after the vectors have been fully populated, with the boundary type of 
     * node being the only thing that is to be identified.
     * 
     * @param lattice_dimensions Extension of each lattice dimension
     * @param nodes Vector of nodes to be checked
     */
    template<>
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

    /**
     * @brief Identifies the propagation direction of each obstacle node based in the general d-dimensional computational domain. Although being
     * a general method, it is not implemented for the general case, as the propagation direction strictly depends on the chosen lattice configuration.
     * @tparam d dimensions of the computational domain
     * @param lattice_dimensions extension of each dimension of the lattice
     * @param nodes obstacle nodes
     */
    template<std::size_t d>
    void identify_obstacle_propagation(const std::array<Eigen::Index, d>& lattice_dimensions,
        std::vector<ObstaclePoint<d>>& nodes,
        const std::vector<Point<d>>& fluids,
        const std::vector<BoundaryPoint<d>>& inlets,
        const std::vector<BoundaryPoint<d>>& outlets)
    {
        return;
    }

    /**
     * @brief Identifies the propagation direction of each obstacle node based on the 2D computational domain. Although being computationally 
     * expensive, the method is used to identify the propagation direction of each obstacle node in the lattice and is needed only once in the generation
     * of lattice, when the directions are identified. Each obstacle node is checked for the presence of fluid nodes in its vicinity, and the propagation direction is set accordingly 
     * by acting on the bitfield of each node.
     * @param lattice_dimensions extension of each dimension of the lattice
     * @param nodes obstacle nodes
     */
    template<>
    void identify_obstacle_propagation(const std::array<Eigen::Index, 2>& lattice_dimensions,
        std::vector<ObstaclePoint<2>>& nodes,
        const std::vector<Point<2>>& fluids,
        const std::vector<BoundaryPoint<2>>& inlets,
        const std::vector<BoundaryPoint<2>>& outlets)
    {

        std::map<std::size_t, std::size_t> d2q9_to_linearized_map{
            {0, 6},
            {1, 2},
            {2, 5},
            {3, 3},
            {4, 0},
            {5, 1},
            {6, 7},
            {7, 4},
            {8, 8}
        };

        for (auto& node : nodes)
        {
            // We know that an obstacle cannot be set at the boundary of the lattice, so we can safely assume that the node is not at the boundary.
            // Thus, we can check the presence of fluid nodes in the direct vicinity of the obstacle with a for loop.
            for (int i = -1; i <= 1; i++)
            {
                for (int j = -1; j <= 1; j++)
                {
                    if (i == 0 && j == 0)
                    {
                        node.directions.set(0, true);
                        continue;
                    }
                    else
                    {
                        // Now if we find that the coordinates of the fluid node match a node that is either in the fluid, inlets, or outlets, we can set the 
                        // corresponding propagation direction to true.
                        const Eigen::Index i_node = node.coords[0] + i;
                        const Eigen::Index j_node = node.coords[1] + j;

                        Point<2> reference_node = {i_node, j_node};
                        const auto fluid_it = std::find(fluids.begin(), fluids.end(), reference_node);
                        std::size_t idx = d2q9_to_linearized_map.at((i+1)*3 + (j+1));
                        if (fluid_it != fluids.end())
                        {
                            node.directions.set(idx, true);
                        }
                        else
                        {
                            const auto inlet_it = std::find(inlets.begin(), inlets.end(), reference_node);
                            if (inlet_it != inlets.end())
                            {
                                node.directions.set(idx, true);
                            }
                            else
                            {
                                const auto outlet_it = std::find(outlets.begin(), outlets.end(), reference_node);
                                if (outlet_it != outlets.end())
                                {
                                    node.directions.set(idx, true);
                                }
                            }
                        }
                    }
                }
            }


        }
    }
}


#endif