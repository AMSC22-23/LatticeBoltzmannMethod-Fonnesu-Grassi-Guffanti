/**
 * @file GenerationTest.cpp
 * @brief Generation test uses the lattice generation infrastructure to build 
 * different possible lattices.
 * 
 * 1. Lid driven cavity lattice.
 */


#include <iostream>

#include "llalbm.hpp"

int main()
{
    using namespace llalbm::util::generation;
    using namespace llalbm::core;
    using namespace llalbm::core::equilibrium;

    ConstructionInfo<2> construction_info;
    
    // Build a domain of 10x10 elements
    construction_info.attach_domain_dimensions({10, 10});
    // Introduce inlets in the first row.
    construction_info.add_nodes_interval({0, 1}, {0, 8}, NonFluidNodeType::INLET);

    auto nodes = construction_info.get_inlet_nodes();
    // print the nodes
    for (auto node : nodes)
    {
        std::cout << "Inlet node: ";
        for (auto coord : node.coords)
        {
            std::cout << coord << " ";
        }
        std::cout << std::endl;
    }
    // Introduce boundaries all around.
    construction_info.add_nodes_interval({0, 0}, {9, 0}, NonFluidNodeType::BOUNDARY);
    construction_info.add_nodes_interval({9, 0}, {9, 9}, NonFluidNodeType::BOUNDARY);
    construction_info.add_nodes_interval({0, 9}, {9, 9}, NonFluidNodeType::BOUNDARY);

    nodes = construction_info.get_boundary_nodes();
    // Introduce boundaries all around.


    using Config = LatticeConfiguration<
        2,
        collisions::BGKCollisionPolicy<2>,
        boundaries::BounceBackPolicy<2>,
        boundaries::BounceBackPolicy<2>,
        boundaries::ZouHePolicy<2>,
        boundaries::ZouHePolicy<2>, 
        initializers::VelocityInitializer<2>,
        equilibrium::DefaultEquilibrium<2>
    >;   

    using Parallel = SerialPolicy<2, Config>;
    
    Lattice<Config, Parallel> Lid;
    std::cout << "Building Lid driven cavity lattice from information" << std::endl;
    build_lattice<2, Config, Parallel>(Lid, 9, construction_info);

    std::ofstream out("file.txt");

    Lid.print_lattice_structure(out, true);

    ConstructionInfo<2> construction_info2_expected_fail;
    build_lattice<2, Config, Parallel>(Lid, 9, construction_info2_expected_fail);

    return 0;
}