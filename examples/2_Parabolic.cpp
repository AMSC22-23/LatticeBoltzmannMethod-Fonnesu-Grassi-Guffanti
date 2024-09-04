#define LLALBM_DEBUG

#include <llalbm.hpp>
#include <iostream>


/**
 * @brief In this second example we are going to show how to simulate a parabolic flow. Here we 
 * want to simulate horizontal flow in a channel with an imposed a parabolic profile. An outlet is on the right and
 * imposes the same profile to the outgoing flow.
 * 
 * In order to do this we can simply modify the function passed as parameter to the initializer. 
 */
int main()
{
    using namespace llalbm::core;
    using namespace llalbm::util::generation;

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
    info.attach_domain_dimensions({50, 200});
    // Then we will place boundaries all along the perimeter
    info.add_perimeter_nodes(NonFluidNodeType::BOUNDARY);
    // We will place an inlet on the left
    info.add_nodes_interval({1, 0}, {48, 0}, NonFluidNodeType::INLET);
    // And an outlet on the right
    info.add_nodes_interval({1, 199}, {48, 199}, NonFluidNodeType::OUTLET);

    // The lattice is ready to be built, and to do that we just need to call build_lattice!
    build_lattice<2, Parallelization>(lattice, 9, info);

    // We can make sure that the lattice has been built correctly by printing the structure to a file
    std::ofstream out("PouseilleStructure.txt");
    lattice.print_lattice_structure(out);

    // Setting up the collider and initializer
    collisions::OMPTRTCollisionPolicy<2>::initialize(0.58, 1.72, 1.0/std::sqrt(3.0));

    // Now we define the inlet velocities
    std::array<std::function<double(double, BoundaryPoint<2>)>, 2> inlets = {
        [](double t, BoundaryPoint<2> p) { 
            // Parabolic flow profile
            const double y = p.coords[0];
            return 0.2*(-(y-25.0)*(y-25.0)+25*25)/(25*25) * (1-std::exp(-t * (0.001)));
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

    // We are done! Let's start the simulation!
    lattice.perform_lbm(3000, 1, 30);
}