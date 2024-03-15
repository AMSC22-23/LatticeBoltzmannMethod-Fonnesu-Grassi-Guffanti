/**
 * @file BounceBackPolicy.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-02-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <vector>
#include <map>
#include <array>
// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "../../utils/loggers/Logger.hpp"
#include "../../utils/aliases.hpp" 

// =======================================


namespace llalbm::core::boundaries
{
    
    using namespace llalbm::util::logger;  
    //using namespace Eigen;  

    template<
        std::size_t dim>
    class BounceBackPolicy
    {
    private:
        /* data */
        Logger log;
    public:
        BounceBackPolicy(/* args */)
        : log("Generic BB" , std::cout)
        {
            log.error("Generic BB not implemented");
            assert(false);
        }
    };

    template<>
    class BounceBackPolicy<2>{
        
        private:
            //std::array<Eigen::Index, 2> lattice_nodes; 
            //std::vector<Point<2>> boundary_nodes;
            Eigen::Index i, j;
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            
        public: 
            /*void load_nodes(std::array<Eigen::Index, 2> &l, std::vector<Point<2>> &b){
                lattice_nodes = l;
                boundary_nodes = b;
            }*/
            void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord)
            {
                auto n_rows = populations.dimensions()[0];
                auto n_cols = populations.dimensions()[1];
                for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++) // * per castare il pointer????
                {
                    i = boundary_coord[bnode].coords[0];
                    j = boundary_coord[bnode].coords[1];

                    p0 = populations(i, j, 0);
                    p1 = populations(i, j, 1);
                    p2 = populations(i, j, 2);
                    p3 = populations(i, j, 3);
                    p4 = populations(i, j, 4);
                    p5 = populations(i, j, 5);
                    p6 = populations(i, j, 6);
                    p7 = populations(i, j, 7);
                    p8 = populations(i, j, 8);

                    populations(i,j,0) = p0;
                    if(j != 0)
                    {
                        populations(i, j-1, 3) = p1;
                    }
                    if(i != n_rows-1)
                    {
                        populations(i+1, j, 4) = p2;
                    }
                    if(j != n_cols-1)
                    {
                        populations(i, j+1, 1) = p3;
                    }
                    if(i != 0)
                    {
                        populations(i-1, j, 2) = p4;
                    }
                    if(i != n_rows-1 && j != 0)
                    {
                        populations(i+1, j-1, 7) = p5;
                    }
                    if(i != n_rows-1 && j != n_cols-1)
                    {
                        populations(i+1, j+1, 8) = p6;
                    }
                    if(i != 0 && j != n_cols-1)
                    {
                        populations(i-1, j+1, 5) = p7;
                    }
                    if(i != 0 && j != 0)
                    {
                        populations(i-1, j-1, 6) = p8;
                    }
                }
            }
    };



} // namespace llalbm::core::boundaries
