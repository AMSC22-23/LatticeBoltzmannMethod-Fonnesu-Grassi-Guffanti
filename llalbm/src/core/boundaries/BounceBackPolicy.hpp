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

                    // populations(i,j,0) = p0;

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

            static void update_boundaries(Tensor<double, 3> &populations, std::vector<ObstaclePoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                const std::size_t obstacles = boundary_coord.size();
                for (std::size_t o_node = 0; o_node < obstacles; ++o_node)
                {
                    const Eigen::Index i = boundary_coord[o_node].coords[0];
                    const Eigen::Index j = boundary_coord[o_node].coords[1];

                    const std::bitset<9> bb_mask = boundary_coord[o_node].directions;

                    if (bb_mask.test(1))
                    {
                        populations(i, j+1, 1) = populations(i, j, 3);
                    }
                    if (bb_mask.test(2))
                    {
                        populations(i-1, j, 2) = populations(i, j, 4);
                    }
                    if (bb_mask.test(3))
                    {
                        populations(i, j-1, 3) = populations(i, j, 1);
                    }
                    if (bb_mask.test(4))
                    {
                        populations(i+1, j, 4) = populations(i, j, 2);
                    }
                    if (bb_mask.test(5))
                    {
                        populations(i-1, j+1, 5) = populations(i, j, 7);
                    }
                    if (bb_mask.test(6))
                    {
                        populations(i-1, j-1, 6) = populations(i, j, 8);
                    }
                    if (bb_mask.test(7))
                    {
                        populations(i+1, j-1, 7) = populations(i, j, 5);
                    }
                    if (bb_mask.test(8))
                    {
                        populations(i+1, j+1, 8) = populations(i, j, 6);
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
                
                #pragma omp parallel for private(i,j,p0,p1,p2,p3,p4,p5,p6,p7,p8)
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

            static void update_boundaries(Tensor<double, 3> &populations, std::vector<ObstaclePoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                const std::size_t obstacles = boundary_coord.size();
                
                #pragma omp parallel for firstprivate(obstacles)
                for (std::size_t o_node = 0; o_node < obstacles; ++o_node)
                {
                    const Eigen::Index i = boundary_coord[o_node].coords[0];
                    const Eigen::Index j = boundary_coord[o_node].coords[1];

                    const std::bitset<9> bb_mask = boundary_coord[o_node].directions;

                    if (bb_mask.test(1))
                    {
                        populations(i, j+1, 1) = populations(i, j, 3);
                    }
                    if (bb_mask.test(2))
                    {
                        populations(i-1, j, 2) = populations(i, j, 4);
                    }
                    if (bb_mask.test(3))
                    {
                        populations(i, j-1, 3) = populations(i, j, 1);
                    }
                    if (bb_mask.test(4))
                    {
                        populations(i+1, j, 4) = populations(i, j, 2);
                    }
                    if (bb_mask.test(5))
                    {
                        populations(i-1, j+1, 5) = populations(i, j, 7);
                    }
                    if (bb_mask.test(6))
                    {
                        populations(i-1, j-1, 6) = populations(i, j, 8);
                    }
                    if (bb_mask.test(7))
                    {
                        populations(i+1, j-1, 7) = populations(i, j, 5);
                    }
                    if (bb_mask.test(8))
                    {
                        populations(i+1, j+1, 8) = populations(i, j, 6);
                    }
                    
                }
            }
    };


    template<
        std::size_t dim>
    class STDExecBounceBackPolicy : public BoundaryPolicyTag, public STDExecTag{
    private:
        /* data */
        Logger log;
    public:
        STDExecBounceBackPolicy(/* args */)
        : log("Generic BB" , std::cout)
        {
            log.error("Generic BB not implemented");
            assert(false);
        }
    };

    template<>
    class STDExecBounceBackPolicy<2> : public BoundaryPolicyTag, public STDExecTag{
            
        public: 
            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {

                Eigen::Index i, j;
                double p0, p1, p2, p3, p4, p5, p6, p7, p8;

                auto n_rows = populations.dimensions()[0];
                auto n_cols = populations.dimensions()[1];
                //for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++) // * per castare il pointer????
                std::for_each(std::execution::par,boundary_coord.begin(),boundary_coord.end(),[&](const auto bnode)
                {
                    i = bnode.coords[0];
                    j = bnode.coords[1];

                    p0 = populations(i, j, 0);
                    p1 = populations(i, j, 1);
                    p2 = populations(i, j, 2);
                    p3 = populations(i, j, 3);
                    p4 = populations(i, j, 4);
                    p5 = populations(i, j, 5);
                    p6 = populations(i, j, 6);
                    p7 = populations(i, j, 7);
                    p8 = populations(i, j, 8);

                    // populations(i,j,0) = p0;

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
                });
            }

            static void update_boundaries(Tensor<double, 3> &populations, std::vector<ObstaclePoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                //const std::size_t obstacles = boundary_coord.size();
                //for (std::size_t o_node = 0; o_node < obstacles; ++o_node)
                std::for_each(std::execution::par,boundary_coord.begin(),boundary_coord.end(),[&](const auto o_node)
                {
                    const Eigen::Index i = o_node.coords[0];
                    const Eigen::Index j = o_node.coords[1];

                    const std::bitset<9> bb_mask = o_node.directions;

                    if (bb_mask.test(1))
                    {
                        populations(i, j+1, 1) = populations(i, j, 3);
                    }
                    if (bb_mask.test(2))
                    {
                        populations(i-1, j, 2) = populations(i, j, 4);
                    }
                    if (bb_mask.test(3))
                    {
                        populations(i, j-1, 3) = populations(i, j, 1);
                    }
                    if (bb_mask.test(4))
                    {
                        populations(i+1, j, 4) = populations(i, j, 2);
                    }
                    if (bb_mask.test(5))
                    {
                        populations(i-1, j+1, 5) = populations(i, j, 7);
                    }
                    if (bb_mask.test(6))
                    {
                        populations(i-1, j-1, 6) = populations(i, j, 8);
                    }
                    if (bb_mask.test(7))
                    {
                        populations(i+1, j-1, 7) = populations(i, j, 5);
                    }
                    if (bb_mask.test(8))
                    {
                        populations(i+1, j+1, 8) = populations(i, j, 6);
                    }
                    
                });
            }
    };


    template<
        std::size_t dim>
    class OpenACCBounceBackPolicy : public BoundaryPolicyTag, public OpenACCTag{
    private:
        /* data */
        Logger log;
    public:
        OpenACCBounceBackPolicy(/* args */)
        : log("Generic BB" , std::cout)
        {
            log.error("Generic BB not implemented");
            assert(false);
        }
    };

    template<>
    class OpenACCBounceBackPolicy<2> : public BoundaryPolicyTag, public OpenACCTag{
            
        public: 
            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {

                Eigen::Index i, j;
                double p0, p1, p2, p3, p4, p5, p6, p7, p8;

                auto n_rows = populations.dimensions()[0];
                auto n_cols = populations.dimensions()[1];
                double* populations_buffer = new double[n_rows * n_cols * 9];

                // Copia i dati in populations_buffer
                for (int i = 0; i < n_rows; ++i) {
                    for (int j = 0; j < n_cols; ++j) {
                        for (int k = 0; k < 9; ++k) {
                            populations_buffer[i * n_cols * 9 + j * 9 + k] = populations(i, j, k);
                        }
                    }
                }

                #pragma acc data copy(populations_buffer[0:n_rows*n_cols*9])
                #pragma acc cache(populations_buffer[0:n_rows*n_cols*9])
                {
                    #pragma acc parallel loop
                    for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++) {
                        int i = boundary_coord[bnode].coords[0];
                        int j = boundary_coord[bnode].coords[1];

                        double* p = &populations_buffer[i * n_cols * 9 + j * 9];
                        double p0 = p[0];
                        double p1 = p[1];
                        double p2 = p[2];
                        double p3 = p[3];
                        double p4 = p[4];
                        double p5 = p[5];
                        double p6 = p[6];
                        double p7 = p[7];
                        double p8 = p[8];

                        // Esegui le modifiche
                        if (j != 0) populations_buffer[i * n_cols * 9 + (j - 1) * 9 + 3] = p1;
                        if (i != n_rows - 1) populations_buffer[(i + 1) * n_cols * 9 + j * 9 + 4] = p2;
                        if (j != n_cols - 1) populations_buffer[i * n_cols * 9 + (j + 1) * 9 + 1] = p3;
                        if (i != 0) populations_buffer[(i - 1) * n_cols * 9 + j * 9 + 2] = p4;
                        if (i != n_rows - 1 && j != 0) populations_buffer[(i + 1) * n_cols * 9 + (j - 1) * 9 + 7] = p5;
                        if (i != n_rows - 1 && j != n_cols - 1) populations_buffer[(i + 1) * n_cols * 9 + (j + 1) * 9 + 8] = p6;
                        if (i != 0 && j != n_cols - 1) populations_buffer[(i - 1) * n_cols * 9 + (j + 1) * 9 + 5] = p7;
                        if (i != 0 && j != 0) populations_buffer[(i - 1) * n_cols * 9 + (j - 1) * 9 + 6] = p8;
                    }
                }

                // Copia i dati modificati indietro in populations
                for (int i = 0; i < n_rows; ++i) {
                    for (int j = 0; j < n_cols; ++j) {
                        for (int k = 0; k < 9; ++k) {
                            populations(i, j, k) = populations_buffer[i * n_cols * 9 + j * 9 + k];
                        }
                    }
                }

                delete[] populations_buffer;

                
                //#pragma acc data copy() 
/*                #pragma acc parallel loop
                for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++)
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

                    // populations(i,j,0) = p0;

                    // Compute the bounce back populations
                    /*
                                j
                              6 2 5
                           i  3 0 1
                              7 4 8
                    
                   
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
                }*/
            }

            static void update_boundaries(Tensor<double, 3> &populations, std::vector<ObstaclePoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                const std::size_t obstacles = boundary_coord.size();

                // Definisci le dimensioni
                auto n_rows = populations.dimensions()[0];
                auto n_cols = populations.dimensions()[1];
                const std::size_t num_directions = 9;

                // Crea un buffer intermedio
                double* p = new double[n_rows * n_cols * num_directions];

                // Copia i dati da populations al buffer
                for (Eigen::Index i = 0; i < n_rows; ++i) {
                    for (Eigen::Index j = 0; j < n_cols; ++j) {
                        for (std::size_t d = 0; d < num_directions; ++d) {
                            p[i * n_cols * num_directions + j * num_directions + d] = populations(i, j, d);
                        }
                    }
                }

                // Mappa il buffer nella GPU
                #pragma acc data copy(p[0:n_rows * n_cols * num_directions])
                #pragma acc cache(p[0:n_rows * n_cols * num_directions])
                {
                    #pragma acc parallel loop
                    for (std::size_t o_node = 0; o_node < obstacles; ++o_node)
                    {
                        const Eigen::Index i = boundary_coord[o_node].coords[0];
                        const Eigen::Index j = boundary_coord[o_node].coords[1];

                        const std::bitset<9> bb_mask = boundary_coord[o_node].directions;

                        if (bb_mask.test(1)) {
                            p[num_directions + 1] = p[num_directions + 3];
                        }
                        if (bb_mask.test(2)) {
                            p[-n_cols * num_directions + 2] = p[num_directions + 4];
                        }
                        if (bb_mask.test(3)) {
                            p[-num_directions] = p[num_directions + 1];
                        }
                        if (bb_mask.test(4)) {
                            p[num_directions] = p[num_directions + 2];
                        }
                        if (bb_mask.test(5)) {
                            p[-n_cols * num_directions + num_directions + 5] = p[num_directions + 7];
                        }
                        if (bb_mask.test(6)) {
                            p[-n_cols * num_directions - num_directions + 6] = p[num_directions + 8];
                        }
                        if (bb_mask.test(7)) {
                            p[n_cols * num_directions - num_directions + 7] = p[num_directions + 5];
                        }
                        if (bb_mask.test(8)) {
                            p[n_cols * num_directions + num_directions + 8] = p[num_directions + 6];
                        }
                    }
                }

                // Copia i dati modificati indietro in populations
                for (Eigen::Index i = 0; i < n_rows; ++i) {
                    for (Eigen::Index j = 0; j < n_cols; ++j) {
                        for (std::size_t d = 0; d < num_directions; ++d) {
                            populations(i, j, d) = p[i * n_cols * num_directions + j * num_directions + d];
                        }
                    }
                }

                delete[] p;

                /*const std::size_t obstacles = boundary_coord.size();

                #pragma acc parallel loop
                for (std::size_t o_node = 0; o_node < obstacles; ++o_node)
                {
                    const Eigen::Index i = boundary_coord[o_node].coords[0];
                    const Eigen::Index j = boundary_coord[o_node].coords[1];

                    const std::bitset<9> bb_mask = boundary_coord[o_node].directions;

                    if (bb_mask.test(1))
                    {
                        populations(i, j+1, 1) = populations(i, j, 3);
                    }
                    if (bb_mask.test(2))
                    {
                        populations(i-1, j, 2) = populations(i, j, 4);
                    }
                    if (bb_mask.test(3))
                    {
                        populations(i, j-1, 3) = populations(i, j, 1);
                    }
                    if (bb_mask.test(4))
                    {
                        populations(i+1, j, 4) = populations(i, j, 2);
                    }
                    if (bb_mask.test(5))
                    {
                        populations(i-1, j+1, 5) = populations(i, j, 7);
                    }
                    if (bb_mask.test(6))
                    {
                        populations(i-1, j-1, 6) = populations(i, j, 8);
                    }
                    if (bb_mask.test(7))
                    {
                        populations(i+1, j-1, 7) = populations(i, j, 5);
                    }
                    if (bb_mask.test(8))
                    {
                        populations(i+1, j+1, 8) = populations(i, j, 6);
                    }
                    
                }*/
            }
    };




} // namespace llalbm::core::boundaries

#endif // LLALBM_BOUNCEBACKPOLICY_HPP