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

#ifndef LLALBM_ZOUHEPOLICY_HPP
#define LLALBM_ZOUHEPOLICY_HPP

#include <vector>
#include <map>
#include <array>
#include <cassert>
// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "../../utils/loggers/Logger.hpp"
#include "../../utils/aliases.hpp" 
#include "../PolicyTypes.hpp"
// =======================================


namespace llalbm::core::boundaries
{
    using namespace llalbm::util::logger;
    template<std::size_t d>
    using Point = Matrix<std::size_t, d, 1>;

    template<std::size_t dim>
    class ZouHePolicy : public BoundaryPolicyTag, public SequentialTag
    {
    private:
        /* data */
        Logger log;
    public:
        ZouHePolicy(/* args */)
        : log("Generic ZH" , std::cout)
        {
            log.error("Generic ZH not implemented");
            assert(false);
        }
    };
    
    template<>
    class ZouHePolicy<2> : public BoundaryPolicyTag, public SequentialTag
    {       
        public: 
            /*void load_nodes(std::array<Eigen::Index, 2> &l, std::vector<Point<2>> &b){
                lattice_nodes = l;
                boundary_nodes = b;
            }*/
            static constexpr double two_thirds = 2.0/3.0;
            static constexpr double one_sixth = 1.0/6.0;
            static constexpr double one_half = 0.5;

            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                
                Eigen::Index i, j;
                double rho, ru, rv;
                
                for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++)
                {
                    i = boundary_coord[bnode].coords[0];
                    j = boundary_coord[bnode].coords[1];
                    switch (boundary_coord[bnode].type)
                    {
                    case TOP_WALL:     
                        rho = (populations(i,j,0) + populations(i,j,1) + populations(i,j,3) + 2.0 * (populations(i,j,2) + populations(i,j,5) + populations(i,j,6))) / (1.0 + global_u(i,j,1)); 
                        global_rho(i,j) = rho;
                        ru = rho * global_u(i,j,0);
                        rv = rho * global_u(i,j,1);

                        populations(i,j,4) = populations(i,j,2) - two_thirds * rv;
                        populations(i,j,7) = populations(i,j,5) - one_sixth * rv - one_half * ru + one_half * (populations(i,j,1) - populations(i,j,3));
                        populations(i,j,8) = populations(i,j,6) - one_sixth * rv + one_half * ru - one_half * (populations(i,j,1) - populations(i,j,3));
                        
                        break;
                    
                    case BOTTOM_WALL:
                        rho = (populations(i,j,0) + populations(i,j,1) + populations(i,j,3) + 2.0 * (populations(i,j,4) + populations(i,j,7) + populations(i,j,8))) / (1.0 - global_u(i,j,1)); 
                        
                        global_rho(i,j) = rho;
                        ru = rho * global_u(i,j,0);
                        rv = rho * global_u(i,j,1);
                        
                        populations(i,j,2) = populations(i,j,4) + two_thirds * rv;
                        populations(i,j,5) = populations(i,j,7) + one_sixth * rv + one_half * ru - one_half * (populations(i,j,1) - populations(i,j,3));
                        populations(i,j,6) = populations(i,j,8) + one_sixth * rv - one_half * ru + one_half * (populations(i,j,1) - populations(i,j,3));

                        break;
                    case LEFT_WALL:
                        rho = (populations(i,j,0) + populations(i,j,2) + populations(i,j,4) + 2.0 * (populations(i,j,3) + populations(i,j,7) + populations(i,j,6))) / (1.0 - global_u(i,j,0)); 

                        global_rho(i,j) = rho;
                        ru = rho * global_u(i,j,0);
                        rv = rho * global_u(i,j,1);

                        populations(i,j,1) = populations(i,j,3) + two_thirds * ru;
                        populations(i,j,5) = populations(i,j,7) + one_sixth * ru + one_half * rv - one_half * (populations(i,j,2) - populations(i,j,4));
                        populations(i,j,8) = populations(i,j,6) + one_sixth * ru - one_half * rv + one_half * (populations(i,j,2) - populations(i,j,4));

                        break;

                    case RIGHT_WALL:
                        rho = (populations(i,j,0) + populations(i,j,2) + populations(i,j,4) + 2.0 * (populations(i,j,1) + populations(i,j,5) + populations(i,j,8))) / (1.0 + global_u(i,j,0)); 
                        
                        global_rho(i,j) = rho;
                        ru = rho * global_u(i,j,0);
                        rv = rho * global_u(i,j,1);

                        populations(i,j,3) = populations(i,j,1) - two_thirds * ru;
                        populations(i,j,7) = populations(i,j,5) - one_sixth * ru - one_half * rv + one_half * (populations(i,j,2) - populations(i,j,4));
                        populations(i,j,6) = populations(i,j,8) - one_sixth * ru + one_half * rv - one_half * (populations(i,j,2) - populations(i,j,4));

                        break;
                    default:
                        break;
                    }
                }
            }
    };
    

    template<std::size_t dim>
    class OMPZouHePolicy : public BoundaryPolicyTag, public OMPTag
    {
    private:
        /* data */
        Logger log;
    public:
        OMPZouHePolicy(/* args */)
        : log("Generic ZH" , std::cout)
        {
            log.error("Generic ZH not implemented");
            assert(false);
        }
    };
    
    template<>
    class OMPZouHePolicy<2> : public BoundaryPolicyTag, public OMPTag
    {       
        public: 
            static constexpr double two_thirds = 2.0/3.0;
            static constexpr double one_sixth = 1.0/6.0;
            static constexpr double one_half = 0.5;

            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                
                Eigen::Index i, j;
                double rho, ru, rv;
                
                #pragma omp parallel for private(i,j)
                for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++)
                {
                    i = boundary_coord[bnode].coords[0];
                    j = boundary_coord[bnode].coords[1];
                    switch (boundary_coord[bnode].type)
                    {
                    case TOP_WALL:     
                        rho = (populations(i,j,0) + populations(i,j,1) + populations(i,j,3) + 2.0 * (populations(i,j,2) + populations(i,j,5) + populations(i,j,6))) / (1.0 + global_u(i,j,1)); 
                        global_rho(i,j) = rho;
                        ru = rho * global_u(i,j,0);
                        rv = rho * global_u(i,j,1);

                        populations(i,j,4) = populations(i,j,2) - two_thirds * rv;
                        populations(i,j,7) = populations(i,j,5) - one_sixth * rv - one_half * ru + one_half * (populations(i,j,1) - populations(i,j,3));
                        populations(i,j,8) = populations(i,j,6) - one_sixth * rv + one_half * ru - one_half * (populations(i,j,1) - populations(i,j,3));
                        
                        break;
                    
                    case BOTTOM_WALL:
                        rho = (populations(i,j,0) + populations(i,j,1) + populations(i,j,3) + 2.0 * (populations(i,j,4) + populations(i,j,7) + populations(i,j,8))) / (1.0 - global_u(i,j,1)); 
                        
                        global_rho(i,j) = rho;
                        ru = rho * global_u(i,j,0);
                        rv = rho * global_u(i,j,1);
                        
                        populations(i,j,2) = populations(i,j,4) + two_thirds * rv;
                        populations(i,j,5) = populations(i,j,7) + one_sixth * rv + one_half * ru - one_half * (populations(i,j,1) - populations(i,j,3));
                        populations(i,j,6) = populations(i,j,8) + one_sixth * rv - one_half * ru + one_half * (populations(i,j,1) - populations(i,j,3));

                        break;
                    case LEFT_WALL:
                        rho = (populations(i,j,0) + populations(i,j,2) + populations(i,j,4) + 2.0 * (populations(i,j,3) + populations(i,j,7) + populations(i,j,6))) / (1.0 - global_u(i,j,0)); 

                        global_rho(i,j) = rho;
                        ru = rho * global_u(i,j,0);
                        rv = rho * global_u(i,j,1);

                        populations(i,j,1) = populations(i,j,3) + two_thirds * ru;
                        populations(i,j,5) = populations(i,j,7) + one_sixth * ru + one_half * rv - one_half * (populations(i,j,2) - populations(i,j,4));
                        populations(i,j,8) = populations(i,j,6) + one_sixth * ru - one_half * rv + one_half * (populations(i,j,2) - populations(i,j,4));

                        break;

                    case RIGHT_WALL:
                        rho = (populations(i,j,0) + populations(i,j,2) + populations(i,j,4) + 2.0 * (populations(i,j,1) + populations(i,j,5) + populations(i,j,8))) / (1.0 + global_u(i,j,0)); 
                        
                        global_rho(i,j) = rho;
                        ru = rho * global_u(i,j,0);
                        rv = rho * global_u(i,j,1);

                        populations(i,j,3) = populations(i,j,1) - two_thirds * ru;
                        populations(i,j,7) = populations(i,j,5) - one_sixth * ru - one_half * rv + one_half * (populations(i,j,2) - populations(i,j,4));
                        populations(i,j,6) = populations(i,j,8) - one_sixth * ru + one_half * rv - one_half * (populations(i,j,2) - populations(i,j,4));

                        break;
                    default:
                        break;
                    }
                }
            }
    };


    template<std::size_t dim>
    class CUDAZouHePolicy : public BoundaryPolicyTag, public CUDATag
    {
    private:
        /* data */
        Logger log;
    public:
        CUDAZouHePolicy(/* args */)
        : log("Generic ZH" , std::cout)
        {
            log.error("Generic ZH not implemented");
            assert(false);
        }
    };

    template<>
    class CUDAZouHePolicy<2> : public BoundaryPolicyTag, public CUDATag
    {       
        public: 
            static constexpr double two_thirds = 2.0/3.0;
            static constexpr double one_sixth = 1.0/6.0;
            static constexpr double one_half = 0.5;

            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {

                Eigen::Index i, j;
                double rho, ru, rv;

                #pragma omp parallel for private(i,j)
                for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++)
                {
                    i = boundary_coord[bnode].coords[0];
                    j = boundary_coord[bnode].coords[1];
                    switch (boundary_coord[bnode].type)
                    {
                    case TOP_WALL:     
                        rho = (populations(i,j,0) + populations(i,j,1) + populations(i,j,3) + 2.0 * (populations(i,j,2) + populations(i,j,5) + populations(i,j,6))) / (1.0 + global_u(i,j,1)); 
                        global_rho(i,j) = rho;
                        ru = rho * global_u(i,j,0);
                        rv = rho * global_u(i,j,1);

                        populations(i,j,4) = populations(i,j,2) - two_thirds * rv;
                        populations(i,j,7) = populations(i,j,5) - one_sixth * rv - one_half * ru + one_half * (populations(i,j,1) - populations(i,j,3));
                        populations(i,j,8) = populations(i,j,6) - one_sixth * rv + one_half * ru - one_half * (populations(i,j,1) - populations(i,j,3));

                        break;

                    case BOTTOM_WALL:
                        rho = (populations(i,j,0) + populations(i,j,1) + populations(i,j,3) + 2.0 * (populations(i,j,4) + populations(i,j,7) + populations(i,j,8))) / (1.0 - global_u(i,j,1)); 

                        global_rho(i,j) = rho;
                        ru = rho * global_u(i,j,0);
                        rv = rho * global_u(i,j,1);

                        populations(i,j,2) = populations(i,j,4) + two_thirds * rv;
                        populations(i,j,5) = populations(i,j,7) + one_sixth * rv + one_half * ru - one_half * (populations(i,j,1) - populations(i,j,3));
                        populations(i,j,6) = populations(i,j,8) + one_sixth * rv - one_half * ru + one_half * (populations(i,j,1) - populations(i,j,3));

                        break;
                    case LEFT_WALL:
                        rho = (populations(i,j,0) + populations(i,j,2) + populations(i,j,4) + 2.0 * (populations(i,j,3) + populations(i,j,7) + populations(i,j,6))) / (1.0 - global_u(i,j,0)); 

                        global_rho(i,j) = rho;
                        ru = rho * global_u(i,j,0);
                        rv = rho * global_u(i,j,1);

                        populations(i,j,1) = populations(i,j,3) + two_thirds * ru;
                        populations(i,j,5) = populations(i,j,7) + one_sixth * ru + one_half * rv - one_half * (populations(i,j,2) - populations(i,j,4));
                        populations(i,j,8) = populations(i,j,6) + one_sixth * ru - one_half * rv + one_half * (populations(i,j,2) - populations(i,j,4));

                        break;

                    case RIGHT_WALL:
                        rho = (populations(i,j,0) + populations(i,j,2) + populations(i,j,4) + 2.0 * (populations(i,j,1) + populations(i,j,5) + populations(i,j,8))) / (1.0 + global_u(i,j,0)); 

                        global_rho(i,j) = rho;
                        ru = rho * global_u(i,j,0);
                        rv = rho * global_u(i,j,1);

                        populations(i,j,3) = populations(i,j,1) - two_thirds * ru;
                        populations(i,j,7) = populations(i,j,5) - one_sixth * ru - one_half * rv + one_half * (populations(i,j,2) - populations(i,j,4));
                        populations(i,j,6) = populations(i,j,8) - one_sixth * ru + one_half * rv - one_half * (populations(i,j,2) - populations(i,j,4));

                        break;
                    default:
                        break;
                    }
                }
            }
    };

    
} // namespace llalbm::core::boundaries

#endif // LLALBM_ZOUHEPOLICY_HPP