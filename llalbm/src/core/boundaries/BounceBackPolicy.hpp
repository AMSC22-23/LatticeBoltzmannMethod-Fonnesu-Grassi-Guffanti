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
            Eigen::Index x, y;
            // double p0;
            double p1, p2, p3, p4, p5, p6, p7, p8;
            
        public: 
            /*void load_nodes(std::array<Eigen::Index, 2> &l, std::vector<Point<2>> &b){
                lattice_nodes = l;
                boundary_nodes = b;
            }*/
            void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> /*global_rho*/, Tensor<double, 3> /*global_u*/){
                
                for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++) // * per castare il pointer????
                {
                    x = boundary_coord[bnode].coords[0];
                    y = boundary_coord[bnode].coords[1];

                    // p0 = populations(x, y, 0);
                    p1 = populations(x, y, 1);
                    p2 = populations(x, y, 2);
                    p3 = populations(x, y, 3);
                    p4 = populations(x, y, 4);
                    p5 = populations(x, y, 5);
                    p6 = populations(x, y, 6);
                    p7 = populations(x, y, 7);
                    p8 = populations(x, y, 8);

                    // populations(x,y,0) = p0;
                    populations(x,y,1) = p3;
                    populations(x,y,2) = p4;
                    populations(x,y,3) = p1;
                    populations(x,y,4) = p2;
                    populations(x,y,5) = p7;
                    populations(x,y,6) = p8;
                    populations(x,y,7) = p5;
                    populations(x,y,8) = p6;

                    /*
                    switch (std::get<1>(boundary_coord[bnode]))
                    {
                    case TOP_WALL:                        
                        populations(x,y,4) = populations(x,y,2);
                        populations(x,y,8) = populations(x,y,6);
                        populations(x,y,7) = populations(x,y,5);
                        break;
                    
                    case LEFT_WALL:
                        populations(x,y,1) = populations(x,y,3);
                        populations(x,y,5) = populations(x,y,7);
                        populations(x,y,8) = populations(x,y,6);
                        break;

                    case RIGHT_WALL:
                        populations(x,y,3) = populations(x,y,1);
                        populations(x,y,7) = populations(x,y,5);
                        populations(x,y,6) = populations(x,y,8);
                        break;

                    case BOTTOM_WALL:
                        populations(x,y,2) = populations(x,y,4);
                        populations(x,y,6) = populations(x,y,8);
                        populations(x,y,5) = populations(x,y,7);
                        break;

                    case RIGHT_BOTTOM_C:
                        populations(x,y,3) = populations(x,y,1);
                        populations(x,y,6) = populations(x,y,8);
                        populations(x,y,2) = populations(x,y,4);
                        break;

                    case LEFT_BOTTOM_C:
                        populations(x,y,2) = populations(x,y,4);
                        populations(x,y,1) = populations(x,y,3);
                        populations(x,y,5) = populations(x,y,7);
                        break;

                    case RIGHT_TOP_C:
                        populations(x,y,3) = populations(x,y,1);
                        populations(x,y,7) = populations(x,y,5);
                        populations(x,y,4) = populations(x,y,2);
                        break;

                    case LEFT_TOP_C:
                        populations(x,y,1) = populations(x,y,3);
                        populations(x,y,8) = populations(x,y,6);
                        populations(x,y,4) = populations(x,y,2);
                        break;
                    }*/
                }
                
            }
    };



} // namespace llalbm::core::boundaries
