#define LLALBM_DEBUG

#include <llalbm.hpp>
#include <iostream>

/**
 * @brief Lid driven cavity test made using the lattice construction infrastructure.
 * 
 * This first example is introductory: it shows how to build a simple lattice using the Lattice Construction Infrastructure
 * provided by LLALBM, and it shows how to perform a simple LBM simulation on it. 
 * 
 * The lattice is a 100x100 lattice with an inlet on the top. The inlet is initialized with horizontal velocity of that increases
 * from 0 to 0.1.
 * 
 * The structure of the lattice is printed to a file called "LidDrivenStructure.txt" and the simulation is performed for 
 * 10000 time steps, with a save every 100 time steps. In this simulation we use the TRT collision model with tau = 0.52,
 * forcing the computation of the second relaxation time to have stability, as given by the magic value.
 * 
 */
int main()
{
    // Core namespace contains the main collision, boundary, initialization and equilibrium policies as well as the lattice
    using namespace llalbm::core;
    // Util namespace contains the lattice generation infrastructure
    using namespace llalbm::util::generation;

    // The first thing to do is to define the Lattice configuration: LLALBM is templated and is built around the concept
    // of policy based design.
    // The Lattice object that we are going to build requires a ParallelizationPolicy template parameter, which requires a 
    // LatticeConfiguration template parameter. LatticeConfiguration contains all the policies that are going to be used to 
    // compute collisions, boundaries, initialization and equilibrium.

    using Configuration = LatticeConfiguration<
        2, // We are using a 2D Lattice
        collisions::OMPTRTCollisionPolicy<2>,
        boundaries::OMPBounceBackPolicy<2>,
        boundaries::OMPBounceBackPolicy<2>,
        boundaries::OMPZouHePolicy<2>,
        boundaries::OMPZouHePolicy<2>,
        initializers::OMPVelocityInitializer<2>,
        equilibrium::OMPDefaultEquilibrium<2>
    >;
    // ... We are using OMP parallelization!

    // Now we need define the type alias for the Parallelization Policy. As we are using OMP, we use the OMPPolicy.
    using Parallelization = OMPPolicy<2, Configuration>;

    // No more type aliases, and now we can build the lattice. As we are going to use the Lattice Construction Infrastructure,
    // we build the lattice object without any parameters. The lattice object will be then reinitialized with all the necessary
    // construction information. 
    Lattice<Parallelization> lattice;

    // Now we need to define the construction information, with the ConstructionInfo object. For now we will not be focusing
    // on the generation of obstacles, that will be analyzed in the next examples.
    ConstructionInfo<2> info;

    // First we need to define the computational domain
    info.attach_domain_dimensions({129, 129});
    // Then we will place boundaries all along the perimeter
    info.add_perimeter_nodes(NonFluidNodeType::BOUNDARY);
    // We will place an inlet on the top
    info.add_nodes_interval({0, 1}, {0, 127}, NonFluidNodeType::INLET);

    // The lattice is ready to be built, and to do that we just need to call build_lattice!
    build_lattice<2, Parallelization>(lattice, 9, info);

    // We can make sure that the lattice has been built correctly by printing the structure to a file
    std::ofstream out("LidDrivenStructure.txt");
    lattice.print_lattice_structure(out);

    // At this point the lattice is ready to be used. We just need to define the simulation parameters.

    // We start by initializing the collider
    collisions::OMPTRTCollisionPolicy<2>::initialize(0.58, 1.0, 1.0/std::sqrt(3.0));
    // But I am not sure that in this case we will have stability... Let's force the computation of the second relaxation time
    collisions::OMPTRTCollisionPolicy<2>::enforce_magic_parameter();

    // Now we define the inlet velocities
    std::array<std::function<double(double, BoundaryPoint<2>)>, 2> inlets = {
        [](double t, BoundaryPoint<2> p) { return 0.15*(1.0 - std::exp(-t * (0.001))); },
        [](double t, BoundaryPoint<2> p) { return 0.0; }
    };

    std::array<std::function<double(double, BoundaryPoint<2>)>, 2> outlets = {
        [](double t, BoundaryPoint<2> p) { return 0.0; },
        [](double t, BoundaryPoint<2> p) { return 0.0; }
    };
    
    initializers::OMPVelocityInitializer<2>::attach_update_functions(inlets, outlets);

    // We are done! Let's start the simulation!
    lattice.perform_lbm(35000, 1, 100);
}