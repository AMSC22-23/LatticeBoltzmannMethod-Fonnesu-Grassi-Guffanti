#define LLALBM_DEBUG

#include <llalbm.hpp>
#include <iostream>


/**
 * @brief In this seventh example we will simulate the flow around a flat-bottom wing profile.
 * To analyze the results we will be instantiating a FlowAnalyzer object and attaching it to the lattice to see how 
 * changes in flow intensity affect lift and drag.
 * 
 */
int main()
{
    using namespace llalbm::core;
    using namespace llalbm::util::generation;
    using namespace llalbm::analysis;

    using Configuration = LatticeConfiguration<
        2, // We are using a 2D Lattice
        collisions::OMPTRTCollisionPolicy<2>,
        boundaries::OMPBounceBackPolicy<2>,
        boundaries::OMPBounceBackPolicy<2>,
        boundaries::OMPZouHePolicy<2>,
        boundaries::OMPOpenBoundaryPolicy<2>,
        initializers::OMPVelocityInitializer<2>,
        equilibrium::OMPDefaultEquilibrium<2>
    >;

    using Parallelization = OMPPolicy<2, Configuration>;

    Lattice<Parallelization> lattice;

    ConstructionInfo<2> info;
    info.attach_domain_dimensions({100, 400});
    info.add_perimeter_nodes(NonFluidNodeType::BOUNDARY);
    info.add_nodes_interval({1, 0}, {98, 0}, NonFluidNodeType::INLET);
    info.add_nodes_interval({1, 399}, {98, 399}, NonFluidNodeType::OUTLET);
    // We will read the wing profile from a file that contains the coordinates of the profile
    info.read_obstacle_from_file("../examples/wing1.txt");
    build_lattice<2, Parallelization>(lattice, 9, info);

    // We can make sure that the lattice has been built correctly by printing the structure to a file
    std::ofstream out("WingStructure1.txt");
    lattice.print_lattice_structure(out);

    // Setting up the collider and initializer
    {
        collisions::OMPTRTCollisionPolicy<2>::initialize(0.6, 1.67, 1.0/std::sqrt(3.0));

        std::array<std::function<double(double, BoundaryPoint<2>)>, 2> inlets = {
            [](double t, BoundaryPoint<2> p) { 
                // Parabolic flow profile
                const double y = p.coords[0];
                return 0.15 * (1-std::exp(-t * (0.001)));
            },
            [](double t, BoundaryPoint<2> p) { return 0.0; }
        };

        std::array<std::function<double(double, BoundaryPoint<2>)>, 2> outlets = {
            [](double t, BoundaryPoint<2> p) { 
                return 0.0;
            },
            [](double t, BoundaryPoint<2> p) { return 0.0; }
        };

        initializers::OMPVelocityInitializer<2>::attach_update_functions(inlets, outlets);
    }

    // We analyze the flow during the simulation, as described in 6_FlowAnalysis.cpp

    std::shared_ptr<FlowAnalyzer<2>> analysis = std::make_shared<FlowAnalyzer<2>>();
    analysis->add_point_vec(lattice.get_obstacle_nodes());
    analysis->set_should_save(true);
    analysis->set_iterations_between_save(100);
    lattice.attach_flow_analyzer(analysis);
    lattice.perform_lbm(8000, 1, 100);
    return 0;
}