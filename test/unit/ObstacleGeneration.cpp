#include "llalbm.hpp"
#include <iostream>

int main()
{
    using namespace llalbm::core;
    using namespace llalbm::util;

    using Config = LatticeConfiguration<
        2,
        collisions::BGKCollisionPolicy<2>,
        boundaries::BounceBackPolicy<2>,
        boundaries::BounceBackPolicy<2>,
        boundaries::ZouHePolicy<2>,
        boundaries::ZouHePolicy<2>,
        initializers::VelocityInitializer<2>,
        equilibrium::DefaultEquilibrium<2>>;

    using Parallel = SerialPolicy<2, Config>;

    llalbm::core::Lattice<Config, Parallel> lattice;
    generation::ConstructionInfo<2> info;

    info.attach_domain_dimensions({10, 10});
    info.add_perimeter_nodes(generation::NonFluidNodeType::BOUNDARY);
    info.add_obstacle_hyper_rectangle({2,2}, {3,2});
    info.add_obstacle_hyper_square({6,6}, 2);

    generation::build_lattice<2, Config>(lattice, 9, info);

    std::ofstream file("lattice.txt");
    lattice.print_lattice_structure(file);

    return 0;
}