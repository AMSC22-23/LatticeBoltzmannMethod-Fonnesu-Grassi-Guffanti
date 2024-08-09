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

#ifndef LLALBM_BOUNCEBACKPOLICY_HPP
#define LLALBM_BOUNCEBACKPOLICY_HPP

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
#include "../PolicyTypes.hpp"
// =======================================


namespace llalbm::core::boundaries
{
    
    using namespace llalbm::util::logger;
    using namespace llalbm::core;  
    //using namespace Eigen;  

    template<
        std::size_t dim>
    class BounceBackPolicy : public BoundaryPolicyTag, public SequentialTag{
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
    class BounceBackPolicy<2> : public BoundaryPolicyTag, public SequentialTag{
            
        public: 
            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {

                Eigen::Index i, j;
                double p0, p1, p2, p3, p4, p5, p6, p7, p8;

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

                    // Compute the bounce back populations
                    /*
                                j
                              6 2 5
                           i  3 0 1
                              7 4 8
                    */
                   
                    // If the node is not on the left wall we can propagate to the left
                    if(j != 0)
                    {
                        populations(i, j-1, 3) = p1;
                    }
                    // If the node is not on the bottom wall we can propagate downwards
                    if(i != n_rows-1)
                    {
                        populations(i+1, j, 4) = p2;
                    }
                    // If the node is not on the right wall we can propagate to the right
                    if(j != n_cols-1)
                    {
                        populations(i, j+1, 1) = p3;
                    }
                    // If the node is not on the top wall we can propagate upwards
                    if(i != 0)
                    {
                        populations(i-1, j, 2) = p4;
                    }
                    // If the node is not on the bottom left corner we can propagate to the bottom left
                    if(i != n_rows-1 && j != 0)
                    {
                        populations(i+1, j-1, 7) = p5;
                    }
                    // If the node is not on the bottom right corner we can propagate to the bottom right
                    if(i != n_rows-1 && j != n_cols-1)
                    {
                        populations(i+1, j+1, 8) = p6;
                    }
                    // If the node is not on the top right corner we can propagate to the top right
                    if(i != 0 && j != n_cols-1)
                    {
                        populations(i-1, j+1, 5) = p7;
                    }
                    // If the node is not on the top left corner we can propagate to the top left
                    if(i != 0 && j != 0)
                    {
                        populations(i-1, j-1, 6) = p8;
                    }
                }
            }
    };


    template<
        std::size_t dim>
    class OMPBounceBackPolicy : public BoundaryPolicyTag, public OMPTag{
    private:
        /* data */
       Logger log;
    public:
        OMPBounceBackPolicy(/* args */)
        : log("Generic BB" , std::cout)
        {
            log.error("Generic BB not implemented");
            assert(false);
        }
    };

    template<>
    class OMPBounceBackPolicy<2> : public BoundaryPolicyTag, public OMPTag{
            
        public: 
            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {

                Eigen::Index i, j;
                double p0, p1, p2, p3, p4, p5, p6, p7, p8;

                auto n_rows = populations.dimensions()[0];
                auto n_cols = populations.dimensions()[1];
                #pragma omp parallel for 
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

                    // Compute the bounce back populations
                    /*
                                j
                              6 2 5
                           i  3 0 1
                              7 4 8
                    */
                   
                    // If the node is not on the left wall we can propagate to the left
                    if(j != 0)
                    {
                        populations(i, j-1, 3) = p1;
                    }
                    // If the node is not on the bottom wall we can propagate downwards
                    if(i != n_rows-1)
                    {
                        populations(i+1, j, 4) = p2;
                    }
                    // If the node is not on the right wall we can propagate to the right
                    if(j != n_cols-1)
                    {
                        populations(i, j+1, 1) = p3;
                    }
                    // If the node is not on the top wall we can propagate upwards
                    if(i != 0)
                    {
                        populations(i-1, j, 2) = p4;
                    }
                    // If the node is not on the bottom left corner we can propagate to the bottom left
                    if(i != n_rows-1 && j != 0)
                    {
                        populations(i+1, j-1, 7) = p5;
                    }
                    // If the node is not on the bottom right corner we can propagate to the bottom right
                    if(i != n_rows-1 && j != n_cols-1)
                    {
                        populations(i+1, j+1, 8) = p6;
                    }
                    // If the node is not on the top right corner we can propagate to the top right
                    if(i != 0 && j != n_cols-1)
                    {
                        populations(i-1, j+1, 5) = p7;
                    }
                    // If the node is not on the top left corner we can propagate to the top left
                    if(i != 0 && j != 0)
                    {
                        populations(i-1, j-1, 6) = p8;
                    }
                }
            }
    };





} // namespace llalbm::core::boundaries

#endif // LLALBM_BOUNCEBACKPOLICY_HPP