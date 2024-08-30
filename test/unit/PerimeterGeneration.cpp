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

    llalbm::core::Lattice<Parallel> lattice;
    generation::ConstructionInfo<2> info;

    info.attach_domain_dimensions({100, 100});
    info.add_perimeter_nodes(generation::NonFluidNodeType::BOUNDARY);

    generation::build_lattice<2, Parallel>(lattice, 9, info);

    std::ofstream file("lattice.txt");
    lattice.print_lattice_structure(file);

    return 0;
}