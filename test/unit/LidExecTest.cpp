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

    llalbm::core::Lattice<Config, Parallel> Lid("../test/assets/lid-lattice.txt", 9,std::cout);

    
    std::array< std::function<double(double,BoundaryPoint<2>)>,2> VelocityFunctions;
    std::array< std::function<double(double,BoundaryPoint<2>)>,2> Outlets;

    VelocityFunctions[0] = [](double time, BoundaryPoint<2> Point){return 0.2/*Da vedere*/;};
    VelocityFunctions[1] = [](double time, BoundaryPoint<2> Point){return 0.0;};
    
    initializers::VelocityInitializer<2>::attach_update_functions(VelocityFunctions,Outlets);

    collisions::BGKCollisionPolicy<2>::initialize(0.6);

    Lid.perform_lbm(100);

}