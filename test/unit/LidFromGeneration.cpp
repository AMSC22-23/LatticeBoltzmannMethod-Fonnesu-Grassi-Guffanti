
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
    using Config = LatticeConfiguration<2,collisions::BGKCollisionPolicy<2>,
                            boundaries::BounceBackPolicy<2>,
                            boundaries::BounceBackPolicy<2>,
                            boundaries::ZouHePolicy<2>,
                            boundaries::ZouHePolicy<2>, 
                            initializers::VelocityInitializer<2>>;    

    llalbm::core::Lattice<Config> Lid;

    
    std::array< std::function<double(double,BoundaryPoint<2>)>,2> VelocityFunctions;
    std::array< std::function<double(double,BoundaryPoint<2>)>,2> Outlets;

    VelocityFunctions[0] = [](double time, BoundaryPoint<2> Point){return 0.2*(1.0-std::exp(-((500*500*time)/(2*1000*1000))));};
    VelocityFunctions[1] = [](double time, BoundaryPoint<2> Point){return 0.0;};
    
    initializers::VelocityInitializer<2>::attach_update_functions(VelocityFunctions,Outlets);
    collisions::BGKCollisionPolicy<2>::initialize(0.6);

    generation::ConstructionInfo<2> info;

    info.attach_domain_dimensions({100, 100});

    info.add_nodes_interval({0, 1}, {0, 98}, generation::NonFluidNodeType::INLET);
    info.add_nodes_interval({0, 0}, {99, 0}, generation::NonFluidNodeType::BOUNDARY);
    info.add_nodes_interval({99, 0}, {99, 99}, generation::NonFluidNodeType::BOUNDARY);
    info.add_nodes_interval({0, 99}, {99, 99}, generation::NonFluidNodeType::BOUNDARY);

    generation::build_lattice<2, Config>(Lid, 9, info);

    std::ofstream out("file.txt");
    Lid.print_lattice_structure(out, true);

    Lid.perform_lbm(1500, 1, 10);

}