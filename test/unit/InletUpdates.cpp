#include "llalbm.hpp"
#include "iostream"

int main()
{
    using namespace llalbm::core;


    using Config = LatticeConfiguration<
        2,
        collisions::BGKCollisionPolicy<2>,
        boundaries::BounceBackPolicy<2>,
        boundaries::BounceBackPolicy<2>,
        boundaries::ZouHePolicy<2>,
        boundaries::ZouHePolicy<2>, 
        initializers::VelocityInitializer<2>,
        equilibrium::DefaultEquilibrium<2>
    >;   
    
    using Parallel = SerialPolicy<2, Config>;

    Lattice<Config, Parallel> Lid("../test/assets/lid-lattice.txt", 9,std::cout);
    
    std::array<std::function<double(double,BoundaryPoint<2>)>,2> VelocityFunctions;
    std::array<std::function<double(double,BoundaryPoint<2>)>,2> Outlets;

    VelocityFunctions[0] = [](double time, BoundaryPoint<2> Point){return 0.2*(1.0 - std::exp(-time * time / (2.0 * 1000)));};
    VelocityFunctions[1] = [](double time, BoundaryPoint<2> Point){return 0.0;};
    
    initializers::VelocityInitializer<2>::attach_update_functions(VelocityFunctions,Outlets);

    for (double t = 0; t < 2.0; t += 0.05)
    {
        double val = VelocityFunctions[0](t, BoundaryPoint<2>({0,0}));
        initializers::VelocityInitializer<2>::update_nodes(t, Lid.get_global_u(), Lid.get_global_rho());
        // In 2D the lid driven cavity the inlet nodes are on top
        // So we can just look at the top row of the lattice to see if the inlet initialization is actually having an effect

        // => print the first row of the velocity tensor
        auto tensor = Lid.get_global_u();
        std::cout << "At time step " << t << " the inlet velocity is: " << val << std::endl;
        for (size_t i = 1; i < Lid.get_lattice_dimensions()[0]-1; i++)
        {
            assert(tensor(0,i,0) == val && "ERROR: Inlet velocity not updated correctly");
        }
        std::cout << "Checked. " << std::endl; 
    }

    return 0;
}
