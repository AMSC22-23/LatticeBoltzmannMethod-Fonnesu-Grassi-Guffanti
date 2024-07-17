#include "llalbm.hpp"

int main()
{
    using namespace llalbm;

    const std::string path = "../test/assets/lattice.txt";
    std::vector<Point<2>> fluid_nodes;
    std::vector<BoundaryPoint<2>> boundary_coord;
    std::vector<BoundaryPoint<2>> inlet_nodes_coord;
    std::vector<BoundaryPoint<2>> outlet_nodes_coord;
    std::vector<BoundaryPoint<2>> obstacle_nodes;
    std::array<Eigen::Index, 2> lattice_dimensions;


    llalbm::util::reader::read_lattice_file<2>(
        path,
        fluid_nodes, 
        boundary_coord,
        inlet_nodes_coord,
        outlet_nodes_coord,
        obstacle_nodes,
        lattice_dimensions
    );

    using namespace llalbm::core;
    
    using Parallel = SerialPolicy<2, 
        collisions::BGKCollisionPolicy<2>, 
        boundaries::BounceBackPolicy<2>, 
        boundaries::BounceBackPolicy<2>, 
        boundaries::ZouHePolicy<2>, 
        boundaries::ZouHePolicy<2>, 
        initializers::VelocityInitializer<2>, 
        equilibrium::DefaultEquilibrium<2>
        >;

    using Config = LatticeConfiguration<
        2,
        collisions::BGKCollisionPolicy<2>,
        boundaries::BounceBackPolicy<2>,
        boundaries::BounceBackPolicy<2>,
        boundaries::ZouHePolicy<2>,
        boundaries::ZouHePolicy<2>, 
        initializers::VelocityInitializer<2>,
        Parallel
    >;   

    Lattice<Config> Lid("../test/assets/lid-lattice.txt", 9,std::cout);

    std::ofstream out("file.txt");
    Lid.print_lattice_structure(out);

    return 0;
}