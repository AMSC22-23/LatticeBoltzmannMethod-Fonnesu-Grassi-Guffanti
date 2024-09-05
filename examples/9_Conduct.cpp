#define LLALBM_DEBUG

#include <llalbm.hpp>
#include <iostream>


/**
 * @brief In this fifth example we are going to show how to simulate a flow moving inside a computational domain
 * which represents a porous medium. The medium is produced by randomly distrivbuting circular obstacles in the domain,
 * uing the lattice construction functionalities.
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
    info.attach_domain_dimensions({200, 600});
    // Then we will place boundaries all along the perimeter
    info.add_perimeter_nodes(NonFluidNodeType::BOUNDARY);
    // We will place an inlet on the left
    info.add_nodes_interval({1, 0}, {198, 0}, NonFluidNodeType::INLET);
    // And an outlet on the right
    info.add_nodes_interval({1, 599}, {198, 599}, NonFluidNodeType::OUTLET);
    // To add random spheres we just need to call the add_obstacle_random_spheres function.
    // The path assumes that the executable is being called from the build/ folder
    info.read_obstacle_from_file("../examples/duct2.txt");
    // The lattice is ready to be built, and to do that we just need to call build_lattice!
    build_lattice<2, Parallelization>(lattice, 9, info);

    // We can make sure that the lattice has been built correctly by printing the structure to a file
    std::ofstream out("ductStructure.txt");
    lattice.print_lattice_structure(out);

    // Setting up the collider and initializer
    collisions::OMPTRTCollisionPolicy<2>::initialize(0.7, 0, 1.0/std::sqrt(3.0));
    collisions::OMPTRTCollisionPolicy<2>::enforce_magic_parameter();
    

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

    // We are done! Let's start the simulation!
    // prima 4700
    lattice.perform_lbm(10000, 1, 30);
}