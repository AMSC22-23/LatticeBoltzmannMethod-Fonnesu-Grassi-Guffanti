
/**
 * @file SphereGenerationTest.cpp
 * @author Luca Guffanti
 * @brief Tests the correct generation of a spheric obstacle in a lattice.
 * 
 */

#define LLALBM_DEBUG 1
#define TIMER

#include "llalbm.hpp"
#include "iostream"
#include <omp.h>

int main()
{
    using namespace llalbm::core;
    using namespace llalbm::util;
    

    using Config = LatticeConfiguration<
        2,
        collisions::OMPTRTCollisionPolicy<2>,
        boundaries::OMPBounceBackPolicy<2>,
        boundaries::OMPBounceBackPolicy<2>,
        boundaries::OMPZouHePolicy<2>,
        boundaries::OMPZouHePolicy<2>, 
        initializers::OMPVelocityInitializer<2>,
        equilibrium::OMPDefaultEquilibrium<2>
    >;   

    using Parallel = OMPPolicy<2, Config>;

    llalbm::core::Lattice<Parallel> Lid;

    generation::ConstructionInfo<2> info;

    info.attach_domain_dimensions({10, 10});
    info.add_obstacle_hyper_sphere({5, 5}, 2);
    
    generation::build_lattice<2, Parallel>(Lid, 9, info);
    std::ofstream file("SphereGeneration.txt");
    Lid.print_lattice_structure(file);

}