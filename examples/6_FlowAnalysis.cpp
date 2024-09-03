
#define DURING_SIMULATION

#ifdef DURING_SIMULATION
#ifdef OUTSIDE_SIMULATION
#undef OUTSIDE_SIMULATION
#endif
#endif

#ifdef OUTSIDE_SIMULATION
#ifdef DURING_SIMULATION
#undef DURING_SIMULATION
#endif
#endif

#define LLALBM_DEBUG
#include "llalbm.hpp"
#include <iostream>


/**
 * @brief In this sixth example we are going to show the two ways of performing flow analysis on an obstacle, as described in the documentation,
 * The first option is to perform the analysis during the simulation, and the second is to perform the analysis after the simulation.
 * 
 * If you are interested in the first option, you can define DURING_SIMULATION, otherwise you can define OUTSIDE_SIMULATION at line 1 of this file.
 * In this simulation we will be using the same cylinder as in the Vortex example, but we will be placing it in the middle of the domain. Then, we will attach a
 * rectangle to the cylinder.
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

    // First we need to define the computational domain
    info.attach_domain_dimensions({100, 400});
    // Then we will place boundaries all along the perimeter
    info.add_perimeter_nodes(NonFluidNodeType::BOUNDARY);
    // We will place an inlet on the left
    info.add_nodes_interval({1, 0}, {98, 0}, NonFluidNodeType::INLET);
    // And an outlet on the right
    info.add_nodes_interval({1, 399}, {98, 399}, NonFluidNodeType::OUTLET);
    // To add a cylinder we just need to insert an obstacle in the middle of the domain
    info.add_obstacle_hyper_sphere({50, 50}, 15);
    // And to avoid eventual vortices we will add a rectangle in the middle. The rectangle is built by passing the origin (top-letf corner) 
    // and the exteneion along each direction.
    info.add_obstacle_hyper_rectangle({45, 50}, {10, 70});

    // The lattice is ready to be built, and to do that we just need to call build_lattice!
    build_lattice<2, Parallelization>(lattice, 9, info);

    // We can make sure that the lattice has been built correctly by printing the structure to a file   
    std::ofstream out("FlowAnalysisStructure.txt");
    lattice.print_lattice_structure(out);

    // Setting up the collider and initializer
    {
    collisions::OMPTRTCollisionPolicy<2>::initialize(0.53, 1.8, 1.0/std::sqrt(3.0));

    // Now we define the inlet velocities
    std::array<std::function<double(double, BoundaryPoint<2>)>, 2> inlets = {
        [](double t, BoundaryPoint<2> p) { 
            // Parabolic flow profile
            const double y = p.coords[0];
            return 0.2 * (1-std::exp(-t * (0.001)));
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

#ifdef DURING_SIMULATION
    // To make analysis during the simulation we must instantiate a FlowAnalysis object through a std::shared_ptr
    std::shared_ptr<FlowAnalyzer<2>> analysis = std::make_shared<FlowAnalyzer<2>>();
    // Initialize it
    analysis->add_point_vec(lattice.get_obstacle_nodes());
    analysis->set_should_save(true);
    analysis->set_iterations_between_save(100);
    lattice.attach_flow_analyzer(analysis);
    // Data will be automatically saved to files in the flow_analysis folder
#endif
    // We are done! Let's start the simulation!
    lattice.perform_lbm(4000, 1, 30);

#ifdef OUTSIDE_SIMULATION
    // If we want to perform the analysis after the simulation we can do it by calling the analyze_flow method of the lattice
    FlowAnalyzer<2> analysis2;
    analysis2.add_point_vec(lattice.get_obstacle_nodes());
    analysis2.set_should_save(true);
    std::pair<double, double> res = analysis2.compute_flow_properties(lattice.get_populations());
    // Data will be automatically saved to files in the flow_analysis folder
    // And the results (global lift and drag) are in res variable.
#endif
}