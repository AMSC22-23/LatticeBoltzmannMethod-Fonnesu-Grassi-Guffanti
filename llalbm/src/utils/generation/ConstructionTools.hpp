/**
 * @file ConstructionTools.hpp
 * @author Luca Guffanti
 * @brief Implements tools that convert ConstructionInfo into a Lattice object. 
 * 
 */

#ifndef LLALBM_CONSTRUCTIONTOOLS_HPP
#define LLALBM_CONSTRUCTIONTOOLS_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <set>
// ======================================

// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "../loggers/Logger.hpp"
#include "../aliases.hpp" 
#include "./ConstructionInfo.hpp"
#include "../MultiDimensionalLoop.hpp"
#include "../../core/lattices/Lattice.hpp"
#include "../analysis/NodeIdentification.hpp"

// =======================================

namespace llalbm::util::generation
{

using namespace llalbm::core;
using namespace llalbm::util::generation;
using namespace llalbm::util;



/**
 * @brief Reinitializes a lattice by using information contained in a ConstructionInfo object, passed as a parameter.
 * It's important notice that the lattice must already be built (eventually with a default constructor) before calling this function.
 * Template parameters are the same of the lattice
 * 
 * @param lattice Pre-built lattice to be reinitialized
 * @param info ConstructionInfo object from which information is retrieved
 */
template<
    std::size_t dim                 ,   // Spacial dimensions of the simulation.
    typename CollisionPolicy        ,   // Policy managing the interaction between two fluid nodes.
    typename WallPolicy             ,   // Policy managing the interaction between a fluid node and a wall node.
    typename ObstaclePolicy         ,   // Policy managing the interaction between a fluid node and an internal obstacle.
    typename InletPolicy            ,   // Policy managing the interaction between an internal fluid node and an inlet node.
    typename OutletPolicy           ,   // Policy managing the interaction between an internal fluid node and an outlet policy.s
    typename InitializationPolicy       // Policy managing the initialization of the lattice.
>
void build_lattice(Lattice<dim, CollisionPolicy, WallPolicy, ObstaclePolicy, InletPolicy, OutletPolicy, InitializationPolicy>& lattice, const ConstructionInfo<dim>& info)
{
    // Traduce ConstructionInfo sets into vectors
    Logger logger("ConstructionTools", std::cout);
    logger.info("Converting sets to vectors");


    auto boundary_nodes_set = info.get_boundary_nodes();
    auto inlet_nodes_set = info.get_inlet_nodes();
    auto outlet_nodes_set = info.get_outlet_nodes();
    auto obstacle_nodes_set = info.get_obstacle_nodes();

    std::vector<BoundaryPoint<dim>> boundary_nodes;
    std::vector<BoundaryPoint<dim>> inlet_nodes;
    std::vector<BoundaryPoint<dim>> outlet_nodes;
    std::vector<BoundaryPoint<dim>> obstacle_nodes;

    boundary_nodes.resize(boundary_nodes_set.size());
    inlet_nodes.resize(inlet_nodes_set.size());
    outlet_nodes.resize(outlet_nodes_set.size());
    obstacle_nodes.resize(obstacle_nodes_set.size());

    std::copy(boundary_nodes_set.begin(), boundary_nodes_set.end(), boundary_nodes.begin());
    std::copy(inlet_nodes_set.begin(), inlet_nodes_set.end(), inlet_nodes.begin());
    std::copy(outlet_nodes_set.begin(), outlet_nodes_set.end(), outlet_nodes.begin());
    std::copy(obstacle_nodes_set.begin(), obstacle_nodes_set.end(), obstacle_nodes.begin());

    std::vector<Point<dim>> fluid_nodes;

    logger.info("Computing total elements");
    const Eigen::Index total_elems = std::accumulate(info.get_domain_dimensions().begin(), info.get_domain_dimensions().end(), 1, std::multiplies<Eigen::Index>());

    fluid_nodes.resize(total_elems);
    // Build the lattice fluid nodes coordinates by removing the boundary, inlet, outlet and obstacle nodes
    // from the possible coordinates
    Eigen::array<Eigen::Index, dim> grid_positions;
    Eigen::array<Eigen::Index, dim> dimensions = info.get_domain_dimensions();
    llalbm::util::MultiDimensionalLoop<std::vector<Point<dim>>, Point<dim>, dim>::assign_grid_positions(fluid_nodes, dimensions, grid_positions);
    
    // Remove the boundary nodes from the fluid nodes
    logger.info("Removing boundary nodes from fluid nodes");
    for (auto boundary_node : boundary_nodes)
    {
        auto it = std::find(fluid_nodes.begin(), fluid_nodes.end(), boundary_node);
        if (it != fluid_nodes.end())
        {
            fluid_nodes.erase(it);
        }
    }

    // Remove the inlet nodes from the fluid nodes
    logger.info("Removing inlet nodes from fluid nodes");
    for (auto inlet_node : inlet_nodes)
    {
        auto it = std::find(fluid_nodes.begin(), fluid_nodes.end(), inlet_node);
        if (it != fluid_nodes.end())
        {
            fluid_nodes.erase(it);
        }
    }

    logger.info("Removing outlet nodes from fluid nodes");
    for (auto outlet_node : outlet_nodes)
    {
        auto it = std::find(fluid_nodes.begin(), fluid_nodes.end(), outlet_node);
        if (it != fluid_nodes.end())
        {
            fluid_nodes.erase(it);
        }
    }

    logger.info("Removing obstacle nodes from fluid nodes");
    for (auto obstacle_node : obstacle_nodes)
    {
        auto it = std::find(fluid_nodes.begin(), fluid_nodes.end(), obstacle_node);
        if (it != fluid_nodes.end())
        {
            fluid_nodes.erase(it);
        }
    }

    // Last but not least, identify the inlet and outlet nodes
    logger.info("Identifying inlet and outlet nodes");
    identify_node_type(info.get_domain_dimensions(), inlet_nodes);
    identify_node_type(info.get_domain_dimensions(), outlet_nodes);

    // Finally, reinitialize the lattice.
    auto domain = info.get_domain_dimensions();
    lattice.reinit(
        domain,
        fluid_nodes,
        boundary_nodes,
        inlet_nodes,
        outlet_nodes,
        obstacle_nodes
    );
}

};

#endif