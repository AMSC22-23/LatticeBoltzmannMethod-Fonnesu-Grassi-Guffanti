
/**
 * @file LidFromGeneration.cpp
 * @author Luca Guffanti
 * @brief Lid driven cavity test made using the lattice construction infrastructure.
 * 
 */
#include "llalbm.hpp"
#include "iostream"

int main()
{
    using namespace llalbm::core;
    using namespace llalbm::util;
    llalbm::core::Lattice<2,collisions::BGKCollisionPolicy<2>,
                            boundaries::BounceBackPolicy<2>,
                            boundaries::BounceBackPolicy<2>,
                            boundaries::ZouHePolicy<2>,
                            boundaries::ZouHePolicy<2>, 
                            initializers::VelocityInitializer<2>> Lid;

    
    std::array< std::function<double(double,BoundaryPoint<2>)>,2> VelocityFunctions;
    std::array< std::function<double(double,BoundaryPoint<2>)>,2> Outlets;

    VelocityFunctions[0] = [](double time, BoundaryPoint<2> Point){return 0.2/*Da vedere*/;};
    VelocityFunctions[1] = [](double time, BoundaryPoint<2> Point){return 0.0;};
    
    initializers::VelocityInitializer<2>::attach_update_functions(VelocityFunctions,Outlets);
    collisions::BGKCollisionPolicy<2>::initialize(0.6);

    generation::ConstructionInfo<2> info;

    info.attach_domain_dimensions({10, 10});

    info.add_nodes_interval({0, 1}, {0, 8}, generation::NonFluidNodeType::INLET);
    info.add_nodes_interval({0, 0}, {9, 0}, generation::NonFluidNodeType::BOUNDARY);
    info.add_nodes_interval({9, 0}, {9, 9}, generation::NonFluidNodeType::BOUNDARY);
    info.add_nodes_interval({0, 9}, {9, 9}, generation::NonFluidNodeType::BOUNDARY);

    generation::build_lattice<2,collisions::BGKCollisionPolicy<2>,
                                boundaries::BounceBackPolicy<2>,
                                boundaries::BounceBackPolicy<2>,
                                boundaries::ZouHePolicy<2>,
                                boundaries::ZouHePolicy<2>,
                                initializers::VelocityInitializer<2>>(Lid, 9, info);

    Lid.perform_lbm(100);

}