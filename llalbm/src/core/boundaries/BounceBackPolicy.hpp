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
#include "Logger.hpp"
#include "aliases.hpp" 

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
            
        public: 
            /*void load_nodes(std::array<Eigen::Index, 2> &l, std::vector<Point<2>> &b){
                lattice_nodes = l;
                boundary_nodes = b;
            }*/
            void update_boundaries(Tensor<double, 3> &populations, std::vector<boundaryPoint<2>> &boundary_coord){
                for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++) // * per castare il pointer????
                {
                    double x = std::get<0>(boundary_coord[bnode])[0];
                    double y = std::get<0>(boundary_coord[bnode])[1];
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
                    }
                }
                
            }
    };



} // namespace llalbm::core::boundaries
