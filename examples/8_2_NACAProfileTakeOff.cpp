#define LLALBM_DEBUG

#include <llalbm.hpp>
#include <iostream>


/**
 * @brief In this eighth example we will simulate the flow around a NACA profile, in a takeoff context. To simulate takeoff we need
 * to simulate the change in angle of attack of the wing, done by introducing a time-dependent vertical velocity component.
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
    info.attach_domain_dimensions({200, 600});
    info.add_perimeter_nodes(NonFluidNodeType::BOUNDARY);
    info.add_nodes_interval({1, 0}, {198, 0}, NonFluidNodeType::INLET);
    info.add_nodes_interval({1, 599}, {198, 599}, NonFluidNodeType::OUTLET);
    // We will read the wing profile from a file that contains the coordinates of the profile
    info.read_obstacle_from_file("../examples/NACA-foil-takeoff.txt");
    build_lattice<2, Parallelization>(lattice, 9, info);

    // We can make sure that the lattice has been built correctly by printing the structure to a file
    std::ofstream out("NACA-foil-takeoffStructure.txt");
    lattice.print_lattice_structure(out);

    // Setting up the collider and initializer
    {
        collisions::OMPTRTCollisionPolicy<2>::initialize(0.58, 1.4, 1.0/std::sqrt(3.0));

        std::array<std::function<double(double, BoundaryPoint<2>)>, 2> inlets = {
            [](double t, BoundaryPoint<2> p) { 
                // Parabolic flow profile
                const double y = p.coords[0];
                return 0.15 * (1-std::exp(-t * (0.001)));
            },
            [](double t, BoundaryPoint<2> p) { 
                const double y = p.coords[0];
                if (y > 100 && t > 1000) {
                    return 0.15 * (1-std::exp(-t * (0.001)));
                }
                return 0.0;
             }
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