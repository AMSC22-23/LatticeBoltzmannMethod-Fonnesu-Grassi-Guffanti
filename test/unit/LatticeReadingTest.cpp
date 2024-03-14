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

    return 0;
}