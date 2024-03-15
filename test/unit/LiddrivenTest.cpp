#include "llalbm.hpp"
#include "iostream"

int main()
{
    using namespace llalbm::core;
    llalbm::core::Lattice<2,collisions::BGKCollisionPolicy<2>,
                            boundaries::BounceBackPolicy<2>,
                            boundaries::BounceBackPolicy<2>,
                            boundaries::ZouHePolicy<2>,
                            boundaries::ZouHePolicy<2>, initializers::VelocityInitializer<2>> Lid("../test/assets/lid-lattice.txt", 9,std::cout);

    
    std::array< std::function<double(std::size_t,BoundaryPoint<2>)>,2> VelocityFunctions;
    std::array< std::function<double(std::size_t,BoundaryPoint<2>)>,2> Outlets;

    VelocityFunctions[0] = [](std::size_t time, BoundaryPoint<2> Point){return 1/*Da vedere*/;};
    VelocityFunctions[1] = [](std::size_t time, BoundaryPoint<2> Point){return 0;};
    initializers::VelocityInitializer<2>::attach_update_functions(VelocityFunctions,Outlets);
    initializers::VelocityInitializer<2>::attach_tensor(Lid.get_global_u(),Lid.get_lattice_dimensions());
    
    
    


    return 0;
}
