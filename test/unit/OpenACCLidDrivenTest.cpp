
/**
 * @file LidFromGeneration.cpp
 * @author Luca Guffanti
 * @brief Lid driven cavity test made using the lattice construction infrastructure and the new TRT collision policy.
 * 
 */

//#define LLALBM_DEBUG 1
//#define TIMER

#include "llalbm.hpp"
#include "iostream"
#include <openacc.h>

int main()
{
    using namespace llalbm::core;
    using namespace llalbm::util;
    

    using Config = LatticeConfiguration<
        2,
        collisions::OpenACCTRTCollisionPolicy<2>,
        boundaries::OpenACCBounceBackPolicy<2>,
        boundaries::OpenACCBounceBackPolicy<2>,
        boundaries::OpenACCZouHePolicy<2>,
        boundaries::OpenACCZouHePolicy<2>, 
        initializers::OpenACCVelocityInitializer<2>,
        equilibrium::OpenACCDefaultEquilibrium<2>
    >;   

    using Parallel = OpenACCPolicy<2, Config>;

    llalbm::core::Lattice<Parallel> Lid;

    
    std::array< std::function<double(double,BoundaryPoint<2>)>,2> VelocityFunctions;
    std::array< std::function<double(double,BoundaryPoint<2>)>,2> Outlets;

    VelocityFunctions[1] = [](double time, BoundaryPoint<2> Point){
        if (Point.coords[1] == 0)
            return 0.2*(1.0-std::exp(-((500*500*time)/(2*1000*1000))));
        if (Point.coords[1] == 99)
            return -0.2*(1.0-std::exp(-((500*500*time)/(2*1000*1000))));
        return 0.0;
            };
        VelocityFunctions[0] = [](double time, BoundaryPoint<2> Point){
        if (Point.coords[0] == 0)
            return 0.2*(1.0-std::exp(-((500*500*time)/(2*1000*1000))));
        if (Point.coords[0] == 99)
            return -0.2*(1.0-std::exp(-((500*500*time)/(2*1000*1000))));
        return 0.0;
            };
    
    initializers::OpenACCVelocityInitializer<2>::attach_update_functions(VelocityFunctions,Outlets);

    collisions::OpenACCTRTCollisionPolicy<2>::initialize(0.9, 0.01, 1./std::sqrt(3.0));
    collisions::OpenACCTRTCollisionPolicy<2>::compute_magic_parameter();
    collisions::OpenACCTRTCollisionPolicy<2>::enforce_magic_parameter(1.0/4.0);

    generation::ConstructionInfo<2> info;

    info.attach_domain_dimensions({100, 100});

    // Add a right inlet
    info.add_perimeter_nodes(generation::NonFluidNodeType::BOUNDARY);
    info.add_nodes_interval({0,1}, {0,98}, generation::NonFluidNodeType::INLET);

    generation::build_lattice<2, Parallel>(Lid, 9, info);

    std::ofstream out("file.txt");
    Lid.print_lattice_structure(out, true);

    Lid.perform_lbm(1000, 1, 10);

}