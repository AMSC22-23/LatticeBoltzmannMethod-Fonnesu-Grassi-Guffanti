
/**
 * @file LidFromGeneration.cpp
 * @author Luca Guffanti
 * @brief Tube test with an open boundary as an outlet.
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
        boundaries::OMPOpenBoundaryPolicy<2>, 
        initializers::OMPVelocityInitializer<2>,
        equilibrium::OMPDefaultEquilibrium<2>
    >;   

    using Parallel = OMPPolicy<2, Config>;

    llalbm::core::Lattice<Config, Parallel> Lid;

    
    std::array< std::function<double(double,BoundaryPoint<2>)>,2> VelocityFunctions;
    std::array< std::function<double(double,BoundaryPoint<2>)>,2> Outlets;

    VelocityFunctions[0] = [](double time, BoundaryPoint<2> Point){return 0.2*(-1*(Point.coords[0]-10)*(Point.coords[0]-10) + 100)/100.0 *(1-std::exp(-0.001*time));};
    VelocityFunctions[1] = [](double time, BoundaryPoint<2> Point){return 0.0;};
    
    initializers::OMPVelocityInitializer<2>::attach_update_functions(VelocityFunctions,Outlets);

    collisions::OMPTRTCollisionPolicy<2>::initialize(0.9, 0.01, 1./std::sqrt(3.0));
    collisions::OMPTRTCollisionPolicy<2>::compute_magic_parameter();
    collisions::OMPTRTCollisionPolicy<2>::enforce_magic_parameter(1.0/4.0);

    generation::ConstructionInfo<2> info;

    info.attach_domain_dimensions({20, 100});

    // Add a right inlet
    info.add_perimeter_nodes(generation::NonFluidNodeType::BOUNDARY);
    info.add_nodes_interval({1,0}, {18,0}, generation::NonFluidNodeType::INLET);
    info.add_nodes_interval({1,99}, {18,99}, generation::NonFluidNodeType::OUTLET);
    info.add_obstacle_hyper_rectangle({4, 20}, {7, 15});
    info.add_obstacle_hyper_rectangle({11,20}, {2,8});


    generation::build_lattice<2, Config>(Lid, 9, info);
    std::ofstream file("Open.txt");
    Lid.print_lattice_structure(file);

    Lid.perform_lbm(3000, 1, 10);

}