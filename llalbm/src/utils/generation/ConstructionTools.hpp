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
 * @tparam LatticeConfiguration Configuration of the lattice
 * @param lattice Pre-built lattice to be reinitialized
 * @param info ConstructionInfo object from which information is retrieved
 */
template<
    std::size_t dim,
    typename LatticeConfiguration,
    typename Parallelization      
>
void build_lattice(Lattice<LatticeConfiguration, Parallelization>& lattice, const std::size_t q, const ConstructionInfo<dim>& info)
{
    // Traduce ConstructionInfo sets into vectors
    Logger logger("ConstructionTools", std::cout);
    logger.info("Converting sets to vectors");

    if (!info.are_dimensions_provided())
    {
        logger.error("Domain dimensions are not provided.");
        logger.error("Please provide domain dimensions in the ConstructionInfo object");
        logger.error("Exiting...");
        exit(1);
    }

    if (!info.are_boundaries_provided())
    {
        logger.warn("Boundary nodes are not provided.");
    }

    if (!info.are_inlets_provided())
    {
        logger.warn("Inlet nodes are not provided.");
    }

    if (!info.are_outlets_provided())
    {
        logger.warn("Outlet nodes are not provided.");
    }

    if (!info.are_obstacles_provided())
    {
        logger.warn("Obstacle nodes are not provided.");
    }



    auto boundary_nodes_set = info.get_boundary_nodes();
    auto inlet_nodes_set = info.get_inlet_nodes();
    auto outlet_nodes_set = info.get_outlet_nodes();
    auto obstacle_nodes_set = info.get_obstacle_nodes();

    std::vector<BoundaryPoint<dim>> boundary_nodes;
    std::vector<BoundaryPoint<dim>> inlet_nodes;
    std::vector<BoundaryPoint<dim>> outlet_nodes;
    std::vector<BoundaryPoint<dim>> obstacle_nodes;


    logger.info("There are " + std::to_string(boundary_nodes_set.size()) + " boundary nodes");
    logger.info("There are " + std::to_string(inlet_nodes_set.size()) + " inlet nodes");
    logger.info("There are " + std::to_string(outlet_nodes_set.size()) + " outlet nodes");
    logger.info("There are " + std::to_string(obstacle_nodes_set.size()) + " obstacle nodes");
    
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
    
    logger.info("Computing fluid nodes coordinates");
    llalbm::util::MultiDimensionalLoop<std::vector<Point<dim>>, Point<dim>, dim>::assign_grid_positions(fluid_nodes, dimensions, grid_positions);
    

    // Remove any duplicates from fluid nodes to be sure
    std::sort(fluid_nodes.begin(), fluid_nodes.end());
    fluid_nodes.erase(std::unique(fluid_nodes.begin(), fluid_nodes.end()), fluid_nodes.end());

    logger.info("Removing boundary nodes from fluid nodes");
    for (auto boundary_node : boundary_nodes)
    {
        auto it = std::find(fluid_nodes.begin(), fluid_nodes.end(), boundary_node);
        if (it != fluid_nodes.end())
        {
            fluid_nodes.erase(it);
        }
    }

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
            std::cout << std::endl;

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
        q,
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