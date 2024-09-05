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
#include <execution>
#include <openacc.h>
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

    /**
     * @brief Generic Zou HE technique for calculation populations on boundaries
     * 
     * @tparam dim 
     */
    template<std::size_t dim>
    class ZouHePolicy : public BoundaryPolicyTag, public SequentialTag
    {
    private:
        Logger log;
    public:
        ZouHePolicy()
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
    class STDExecZouHePolicy : public BoundaryPolicyTag, public STDExecTag
    {
    private:
        /* data */
        Logger log;
    public:
        STDExecZouHePolicy(/* args */)
        : log("Generic ZH" , std::cout)
        {
            log.error("Generic ZH not implemented");
            assert(false);
        }
    };
    
    template<>
    class STDExecZouHePolicy<2> : public BoundaryPolicyTag, public STDExecTag
    {       
        public: 
            static constexpr double two_thirds = 2.0/3.0;
            static constexpr double one_sixth = 1.0/6.0;
            static constexpr double one_half = 0.5;

            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                
                Eigen::Index i, j;
                double rho, ru, rv;
                
                std::for_each(std::execution::par,boundary_coord.begin(),boundary_coord.end(),[&](const auto bnode)
                {
                    i = bnode.coords[0];
                    j = bnode.coords[1];
                    switch (bnode.type)
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
                });
            }
    };
    
    template<std::size_t dim>
    class OpenACCZouHePolicy : public BoundaryPolicyTag, public OpenACCTag
    {
    private:
        /* data */
        Logger log;
    public:
        OpenACCZouHePolicy(/* args */)
        : log("Generic ZH" , std::cout)
        {
            log.error("Generic ZH not implemented");
            assert(false);
        }
    };
    
    template<>
    class OpenACCZouHePolicy<2> : public BoundaryPolicyTag, public OpenACCTag
    {       
        public: 
           
            static constexpr double two_thirds = 2.0/3.0;
            static constexpr double one_sixth = 1.0/6.0;
            static constexpr double one_half = 0.5;

            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                
                auto n_rows = populations.dimensions()[0];
                auto n_cols = populations.dimensions()[1];
                const std::size_t num_directions = 9;

               
                double* populations_buffer = new double[n_rows * n_cols * num_directions];

               
                for (Eigen::Index i = 0; i < n_rows; ++i) {
                    for (Eigen::Index j = 0; j < n_cols; ++j) {
                        for (std::size_t d = 0; d < num_directions; ++d) {
                            populations_buffer[i * n_cols * num_directions + j * num_directions + d] = populations(i, j, d);
                        }
                    }
                }

                
                #pragma acc data copy(populations_buffer[0:n_rows * n_cols * num_directions])
                #pragma acc cache(populations_buffer[0:n_rows * n_cols * num_directions])
                {
                    #pragma acc parallel loop
                    for (std::size_t bnode = 0; bnode < boundary_coord.size(); ++bnode)
                    {
                        Eigen::Index i = boundary_coord[bnode].coords[0];
                        Eigen::Index j = boundary_coord[bnode].coords[1];
                        double rho, ru, rv;

                        std::size_t base_idx = i * n_cols * num_directions + j * num_directions;

                        switch (boundary_coord[bnode].type)
                        {
                            case TOP_WALL:
                                rho = (populations_buffer[base_idx + 0] +
                                    populations_buffer[base_idx + 1] +
                                    populations_buffer[base_idx + 3] +
                                    2.0 * (populations_buffer[base_idx + 2] +
                                            populations_buffer[base_idx + 5] +
                                            populations_buffer[base_idx + 6])) /
                                    (1.0 + global_u(i,j,1));
                                global_rho(i,j) = rho;
                                ru = rho * global_u(i,j,0);
                                rv = rho * global_u(i,j,1);

                                populations_buffer[base_idx + 4] = populations_buffer[base_idx + 2] - two_thirds * rv;
                                populations_buffer[base_idx + 7] = populations_buffer[base_idx + 5] - one_sixth * rv - one_half * ru + one_half * (populations_buffer[base_idx + 1] - populations_buffer[base_idx + 3]);
                                populations_buffer[base_idx + 8] = populations_buffer[base_idx + 6] - one_sixth * rv + one_half * ru - one_half * (populations_buffer[base_idx + 1] - populations_buffer[base_idx + 3]);

                                break;

                            case BOTTOM_WALL:
                                rho = (populations_buffer[base_idx + 0] +
                                    populations_buffer[base_idx + 1] +
                                    populations_buffer[base_idx + 3] +
                                    2.0 * (populations_buffer[base_idx + 4] +
                                            populations_buffer[base_idx + 7] +
                                            populations_buffer[base_idx + 8])) /
                                    (1.0 - global_u(i,j,1));
                                global_rho(i,j) = rho;
                                ru = rho * global_u(i,j,0);
                                rv = rho * global_u(i,j,1);

                                populations_buffer[base_idx + 2] = populations_buffer[base_idx + 4] + two_thirds * rv;
                                populations_buffer[base_idx + 5] = populations_buffer[base_idx + 7] + one_sixth * rv + one_half * ru - one_half * (populations_buffer[base_idx + 1] - populations_buffer[base_idx + 3]);
                                populations_buffer[base_idx + 6] = populations_buffer[base_idx + 8] + one_sixth * rv - one_half * ru + one_half * (populations_buffer[base_idx + 1] - populations_buffer[base_idx + 3]);

                                break;

                            case LEFT_WALL:
                                rho = (populations_buffer[base_idx + 0] +
                                    populations_buffer[base_idx + 2] +
                                    populations_buffer[base_idx + 4] +
                                    2.0 * (populations_buffer[base_idx + 3] +
                                            populations_buffer[base_idx + 7] +
                                            populations_buffer[base_idx + 6])) /
                                    (1.0 - global_u(i,j,0));
                                global_rho(i,j) = rho;
                                ru = rho * global_u(i,j,0);
                                rv = rho * global_u(i,j,1);

                                populations_buffer[base_idx + 1] = populations_buffer[base_idx + 3] + two_thirds * ru;
                                populations_buffer[base_idx + 5] = populations_buffer[base_idx + 7] + one_sixth * ru + one_half * rv - one_half * (populations_buffer[base_idx + 2] - populations_buffer[base_idx + 4]);
                                populations_buffer[base_idx + 8] = populations_buffer[base_idx + 6] + one_sixth * ru - one_half * rv + one_half * (populations_buffer[base_idx + 2] - populations_buffer[base_idx + 4]);

                                break;

                            case RIGHT_WALL:
                                rho = (populations_buffer[base_idx + 0] +
                                    populations_buffer[base_idx + 2] +
                                    populations_buffer[base_idx + 4] +
                                    2.0 * (populations_buffer[base_idx + 1] +
                                            populations_buffer[base_idx + 5] +
                                            populations_buffer[base_idx + 8])) /
                                    (1.0 + global_u(i,j,0));
                                global_rho(i,j) = rho;
                                ru = rho * global_u(i,j,0);
                                rv = rho * global_u(i,j,1);

                                populations_buffer[base_idx + 3] = populations_buffer[base_idx + 1] - two_thirds * ru;
                                populations_buffer[base_idx + 7] = populations_buffer[base_idx + 5] - one_sixth * ru - one_half * rv + one_half * (populations_buffer[base_idx + 2] - populations_buffer[base_idx + 4]);
                                populations_buffer[base_idx + 6] = populations_buffer[base_idx + 8] - one_sixth * ru + one_half * rv - one_half * (populations_buffer[base_idx + 2] - populations_buffer[base_idx + 4]);

                                break;

                            default:
                                break;
                        }
                    }
                }

                for (Eigen::Index i = 0; i < n_rows; ++i) {
                    for (Eigen::Index j = 0; j < n_cols; ++j) {
                        for (std::size_t d = 0; d < num_directions; ++d) {
                            populations(i, j, d) = populations_buffer[i * n_cols * num_directions + j * num_directions + d];
                        }
                    }
                }

                delete[] populations_buffer;


            }
    };
    
} // namespace llalbm::core::boundaries

#endif // LLALBM_ZOUHEPOLICY_HPP