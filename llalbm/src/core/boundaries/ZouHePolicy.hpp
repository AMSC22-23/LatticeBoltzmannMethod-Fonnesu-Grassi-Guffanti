/**
 * @file ZouHePolicy.hpp
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
    template<std::size_t d>
    using Point = Matrix<std::size_t, d, 1>;

    template<
        std::size_t dim>
    class ZouHePolicy
    {
    private:
        /* data */
    public:
        ZouHePolicy(/* args */)
        : log("Generic ZH" , std::cout)
        {
            log.error("Generic ZH not implemented");
            assert(false);
        }
    };
    
    template<>
    class ZouHePolicy<2>
    {
        private:
            //std::array<Eigen::Index, 2> lattice_nodes; 
            //std::vector<Point<2>> boundary_nodes;
            Eigen::Index x, y;
            double rho, ru;
            
        public: 
            /*void load_nodes(std::array<Eigen::Index, 2> &l, std::vector<Point<2>> &b){
                lattice_nodes = l;
                boundary_nodes = b;
            }*/
            static constexpr double two_thirds = 2.0/3.0;
            static constexpr double one_sixth = 1.0/6.0;
            static constexpr double one_half = 0.5;

            void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u){
                for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++) // * per castare il pointer????
                {
                    x = boundary_coord[bnode].coords[0];
                    y = boundary_coord[bnode].coords[1];
                    switch (boundary_coord[bnode].type)
                    {
                    case TOP_WALL:        
                        rho = (populations(x,y,0) + populations(x,y,1) + populations(x,y,3) + 2.0 * (populations(x,y,2) + populations(x,y,5) + populations(x,y,6))) / (1.0 + global_u(x,y,1)); 
                        ru = rho * global_u(x,y,1);
                        populations(x,y,4) = populations(x,y,2) - two_thirds * ru;
                        populations(x,y,7) = populations(x,y,5) - one_sixth * ru + one_half * (populations(x,y,1) - populations(x,y,3));
                        populations(x,y,8) = populations(x,y,6) - one_sixth * ru + one_half * (populations(x,y,3) - populations(x,y,1));
                        break;
                    
                    case LEFT_WALL:
                        rho = (populations(x,y,0) + populations(x,y,2) + populations(x,y,4) + 2.0 * (populations(x,y,3) + populations(x,y,7) + populations(x,y,6))) / (1.0 - global_u(x,y,0)); 
                        ru = rho * global_u(x,y,0);
                        populations(x,y,1) = populations(x,y,3) + two_thirds * ru;
                        populations(x,y,5) = populations(x,y,7) + one_sixth * ru - one_half * (populations(x,y,2) - populations(x,y,4));
                        populations(x,y,8) = populations(x,y,6) + one_sixth * ru - one_half * (populations(x,y,4) - populations(x,y,2));
                        break;

                    case RIGHT_WALL:
                        rho = (populations(x,y,0) + populations(x,y,2) + populations(x,y,4) + 2.0 * (populations(x,y,1) + populations(x,y,5) + populations(x,y,8))) / (1.0 + global_u(x,y,0)); 
                        ru = rho * global_u(x,y,0);
                        populations(x,y,3) = populations(x,y,1) - two_thirds * ru;
                        populations(x,y,7) = populations(x,y,5) - one_sixth * ru + one_half * (populations(x,y,2) - populations(x,y,4));
                        populations(x,y,6) = populations(x,y,8) - one_sixth * ru + one_half * (populations(x,y,4) - populations(x,y,2));
                        break;

                    case BOTTOM_WALL:
                        rho = (populations(x,y,0) + populations(x,y,1) + populations(x,y,3) + 2.0 * (populations(x,y,4) + populations(x,y,7) + populations(x,y,8))) / (1.0 - global_u(x,y,1)); 
                        ru = rho * global_u(x,y,1);
                        populations(x,y,2) = populations(x,y,2) + two_thirds * ru;
                        populations(x,y,5) = populations(x,y,7) + one_sixth * ru - one_half * (populations(x,y,1) - populations(x,y,3));
                        populations(x,y,6) = populations(x,y,8) + one_sixth * ru - one_half * (populations(x,y,3) - populations(x,y,1));
                        break;

                    case RIGHT_BOTTOM_C:
                        populations(x,y,3) = populations(x,y,1) - two_thirds * global_rho(x,y) * global_u(x,y,0);
                        populations(x,y,2) = populations(x,y,4) + two_thirds * global_rho(x,y) * global_u(x,y,1);
                        populations(x,y,6) = populations(x,y,8) - one_sixth * global_rho(x,y) * (global_u(x,y,0) - global_u(x,y,1));

                        populations(x,y,5) = 0.0;
                        populations(x,y,7) = 0.0;

                        populations(x,y,0) = global_rho(x,y) 
                                            - populations(x,y,1) - populations(x,y,2) - populations(x,y,3) - populations(x,y,4)
                                            - populations(x,y,5) - populations(x,y,6) - populations(x,y,7) - populations(x,y,8);
                        
                        break;

                    case LEFT_BOTTOM_C:
                        populations(x,y,1) = populations(x,y,3) - two_thirds * global_rho(x,y) * global_u(x,y,0);
                        populations(x,y,2) = populations(x,y,4) - two_thirds * global_rho(x,y) * global_u(x,y,1);
                        populations(x,y,5) = populations(x,y,7) + one_sixth * global_rho(x,y) * (global_u(x,y,0) + global_u(x,y,1));

                        populations(x,y,6) = 0.0;
                        populations(x,y,8) = 0.0;

                        populations(x,y,0) = global_rho(x,y) 
                                            - populations(x,y,1) - populations(x,y,2) - populations(x,y,3) - populations(x,y,4)
                                            - populations(x,y,5) - populations(x,y,6) - populations(x,y,7) - populations(x,y,8);
                        
                        break;

                    case RIGHT_TOP_C:
                        populations(x,y,3) = populations(x,y,1) - two_thirds * global_rho(x,y) * global_u(x,y,0);
                        populations(x,y,4) = populations(x,y,2) - two_thirds * global_rho(x,y) * global_u(x,y,1);
                        populations(x,y,7) = populations(x,y,5) - one_sixth * global_rho(x,y) * (global_u(x,y,0) + global_u(x,y,1));

                        populations(x,y,6) = 0.0;
                        populations(x,y,8) = 0.0;

                        populations(x,y,0) = global_rho(x,y) 
                                            - populations(x,y,1) - populations(x,y,2) - populations(x,y,3) - populations(x,y,4)
                                            - populations(x,y,5) - populations(x,y,6) - populations(x,y,7) - populations(x,y,8);
                        break;

                    case LEFT_TOP_C:
                        populations(x,y,1) = populations(x,y,3) + two_thirds * global_rho(x,y) * global_u(x,y,0);
                        populations(x,y,4) = populations(x,y,2) - two_thirds * global_rho(x,y) * global_u(x,y,1);
                        populations(x,y,8) = populations(x,y,6) + one_sixth * global_rho(x,y) * (global_u(x,y,0) - global_u(x,y,1));

                        populations(x,y,5) = 0.0;
                        populations(x,y,7) = 0.0;

                        populations(x,y,0) = global_rho(x,y) 
                                            - populations(x,y,1) - populations(x,y,2) - populations(x,y,3) - populations(x,y,4)
                                            - populations(x,y,5) - populations(x,y,6) - populations(x,y,7) - populations(x,y,8);
                        break;
                    }
                }
            }
    };
    
} // namespace llalbm::core::boundaries
