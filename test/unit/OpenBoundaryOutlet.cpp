
/**
 * @file OpenBoundaryOutlet.cpp
 * @author Luca Guffanti
 * @brief Tube test with an open boundary as an outlet.
 * 
 */

#define TIMER
#define LLALBM_DEBUG 1


#include "llalbm.hpp"
#include "iostream"
#include <omp.h>

int main()
{
    using namespace llalbm::core;
    using namespace llalbm::util;
    using namespace llalbm::analysis;
    

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

    
    std::array< std::function<double(double,BoundaryPoint<2>)>,2> VelocityFunctions;
    std::array< std::function<double(double,BoundaryPoint<2>)>,2> Outlets;

    VelocityFunctions[0] = [](double time, BoundaryPoint<2> Point){return 0.04*(1-std::exp(-0.001*time));};
    VelocityFunctions[1] = [](double time, BoundaryPoint<2> Point){return 0.0;};
    
    Outlets[0] = [](double time, BoundaryPoint<2> Point){return 0.0;};
    Outlets[1] = [](double time, BoundaryPoint<2> Point){return 0.0;};

    initializers::OMPVelocityInitializer<2>::attach_update_functions(VelocityFunctions,Outlets);

    collisions::OMPTRTCollisionPolicy<2>::initialize(0.7, 0.01, 1./std::sqrt(3.0));
    collisions::OMPTRTCollisionPolicy<2>::compute_magic_parameter();
    collisions::OMPTRTCollisionPolicy<2>::enforce_magic_parameter(1.0/4.0);

    generation::ConstructionInfo<2> info;

    info.attach_domain_dimensions({100, 100});

    // Add a right inlet
    info.add_perimeter_nodes(generation::NonFluidNodeType::BOUNDARY);
    info.add_nodes_interval({1,0}, {98,0}, generation::NonFluidNodeType::INLET);
    info.add_nodes_interval({1,99}, {98, 99}, generation::NonFluidNodeType::OUTLET);
    info.add_obstacle_hyper_sphere({ 50, 50}, 30);

    generation::build_lattice<2, Parallel>(Lid, 9, info);
    std::ofstream file("Open.txt");
    Lid.print_lattice_structure(file);

    std::shared_ptr<FlowAnalyzer<2>> flow_analyzer = std::make_shared<FlowAnalyzer<2>>();
    flow_analyzer->set_iterations_between_save(10);
    flow_analyzer->add_point_vec(Lid.get_obstacle_nodes());
    
    flow_analyzer->print_considered_points();

    Lid.attach_flow_analyzer(flow_analyzer);

    Lid.perform_lbm(1000, 1, 100);
}