#define LLALBM_DEBUG
#include "llalbm.hpp"
#include <iostream>


int main(){
    using namespace llalbm::core;
    using namespace llalbm::util::reader;
    using namespace llalbm::util::generation;

    std::vector<ObstaclePoint<2>> obstacles;

    
    using Config = LatticeConfiguration<
        2,
        collisions::OMPBGKCollisionPolicy<2>,
        boundaries::OMPBounceBackPolicy<2>,
        boundaries::OMPBounceBackPolicy<2>,
        boundaries::OMPZouHePolicy<2>,
        boundaries::OMPZouHePolicy<2>,
        initializers::OMPVelocityInitializer<2>,
        equilibrium::OMPDefaultEquilibrium<2>
        >;

    using Parallelization = OMPPolicy<2, Config>;

    Lattice<Parallelization> lattice;
    ConstructionInfo<2> info;
    info.attach_domain_dimensions({200, 600});
    info.add_perimeter_nodes(NonFluidNodeType::BOUNDARY);
    info.add_nodes_interval({1,0}, {198, 0}, NonFluidNodeType::INLET);
    info.add_nodes_interval({1,599}, {198, 599}, NonFluidNodeType::OUTLET);
    info.read_obstacle_from_file("black_pixel_coords.txt");

    build_lattice<2, Parallelization>(lattice, 9, info);

    collisions::OMPBGKCollisionPolicy<2>::initialize(0.6, 10);

    std::array<std::function<double(double, BoundaryPoint<2>)>, 2> inlets = {
        [](double t, BoundaryPoint<2> p) { return 0.15*(1.0 - std::exp(-t * (0.001))); },
        [](double t, BoundaryPoint<2> p) { return 0.0; }
    };

    std::array<std::function<double(double, BoundaryPoint<2>)>, 2> outlets = {
        [](double t, BoundaryPoint<2> p) { return 0.0; },
        [](double t, BoundaryPoint<2> p) { return 0.0; }
    };

    initializers::OMPVelocityInitializer<2>::attach_update_functions(inlets, outlets);

    lattice.perform_lbm(400, 1, 100);
    FlowAnalyzer<2> flow_analyzer;
    flow_analyzer.add_point_vec(lattice.get_obstacle_nodes());
    auto data = flow_analyzer.compute_flow_properties(lattice.get_populations());
    std::cout << "Total lift: " << data.first << std::endl;
    std::cout << "Total drag: " << data.second << std::endl;


}