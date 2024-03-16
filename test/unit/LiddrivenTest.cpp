#include "llalbm.hpp"
#include "iostream"

int main()
{
    using namespace llalbm::core;
    llalbm::core::Lattice<2,collisions::BGKCollisionPolicy<2>,
                            boundaries::BounceBackPolicy<2>,
                            boundaries::BounceBackPolicy<2>,
                            boundaries::ZouHePolicy<2>,
                            boundaries::ZouHePolicy<2>, 
                            initializers::VelocityInitializer<2>> Lid("../test/assets/lid-lattice.txt", 9,std::cout);

    
    std::array< std::function<double(std::size_t,BoundaryPoint<2>)>,2> VelocityFunctions;
    std::array< std::function<double(std::size_t,BoundaryPoint<2>)>,2> Outlets;

    VelocityFunctions[0] = [](std::size_t time, BoundaryPoint<2> Point){return 1.0/*Da vedere*/;};
    VelocityFunctions[1] = [](std::size_t time, BoundaryPoint<2> Point){return 0.0;};
    
    initializers::VelocityInitializer<2>::attach_update_functions(VelocityFunctions,Outlets);

    initializers::VelocityInitializer<2>::update_nodes(0, Lid.get_global_u(), Lid.get_global_rho());
    // Get the list of inlets from the lattice
    auto inlets = Lid.get_inlet_nodes();
    // Print the list of inlets
    for (auto inlet : inlets)
    {
        std::cout << "Inlet node: ";
        for (auto coord : inlet.coords)
        {
            std::cout << coord << " ";
        }
        std::cout << std::endl;
    }

    // And do the same for the outlets
    auto outlets = Lid.get_outlet_nodes();
    for (auto outlet : outlets)
    {
        std::cout << "Outlet node: ";
        for (auto coord : outlet.coords)
        {
            std::cout << coord << " ";
        }
        std::cout << std::endl;
    }

    // And for the boundary nodes
    auto boundaries = Lid.get_boundary_nodes();
    for (auto boundary : boundaries)
    {
        std::cout << "Boundary node: ";
        for (auto coord : boundary.coords)
        {
            std::cout << coord << " ";
        }
        std::cout << std::endl;
    }

    // For the obstacle nodes
    auto obstacles = Lid.get_obstacle_nodes();
    for (auto obstacle : obstacles)
    {
        std::cout << "Obstacle node: ";
        for (auto coord : obstacle.coords)
        {
            std::cout << coord << " ";
        }
        std::cout << std::endl;
    }

    // And for the fluids
    auto fluids = Lid.get_fluid_nodes();
    for (auto fluid : fluids)
    {
        std::cout << "Fluid node: ";
        for (auto coord : fluid.coords)
        {
            std::cout << coord << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}
