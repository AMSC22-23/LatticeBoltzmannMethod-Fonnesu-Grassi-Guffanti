#ifndef LLALBM_DEFAULT_EQUILIBRIUM_HPP
#define LLALBM_DEFAULT_EQUILIBRIUM_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
#include <execution>
#include <openacc.h>
// ======================================

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

namespace llalbm::core::equilibrium
{

    template<
        std::size_t dim
        >
    class DefaultEquilibrium : SequentialTag
    {};

    template<>
    class DefaultEquilibrium<2> : SequentialTag {    
    
    public:
        static void calc_equilibrium(std::vector<Point<2>> &fluid_nodes, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &global_u, Tensor<double,2> &global_rho){
            double rho,ux,uy,u2;
            constexpr double two9 = (2.0/9.0);
            constexpr double one18 = (1.0/18.0);
            constexpr double one36 = (1.0/36.0);
            
            for (size_t k = 0; k < fluid_nodes.size(); k++)
            {
                size_t i = fluid_nodes[k].coords[0];
                size_t j = fluid_nodes[k].coords[1];
                rho = global_rho(i,j);
                ux = global_u(i,j,0);
                uy = global_u(i,j,1);
                u2 = ux * ux + uy * uy;
                equilibrium_populations(i,j,0) =  two9 * rho * (2 - 3*u2);
                equilibrium_populations(i,j,1) =  one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                equilibrium_populations(i,j,2) =  one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                equilibrium_populations(i,j,3) =  one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                equilibrium_populations(i,j,4) =  one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                equilibrium_populations(i,j,5) = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                equilibrium_populations(i,j,6) = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                equilibrium_populations(i,j,7) = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                equilibrium_populations(i,j,8) = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
            }            
        }

        static void calc_equilibrium(std::vector<BoundaryPoint<2>> &fluid_nodes, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &global_u, Tensor<double,2> &global_rho){
            double rho,ux,uy,u2;
            constexpr double two9 = (2.0/9.0);
            constexpr double one18 = (1.0/18.0);
            constexpr double one36 = (1.0/36.0);
            
            for (size_t k = 0; k < fluid_nodes.size(); k++)
            {
                size_t i = fluid_nodes[k].coords[0];
                size_t j = fluid_nodes[k].coords[1];
                rho = global_rho(i,j);
                ux = global_u(i,j,0);
                uy = global_u(i,j,1);
                u2 = ux * ux + uy * uy;
                equilibrium_populations(i,j,0) =  two9 * rho * (2 - 3*u2);
                equilibrium_populations(i,j,1) =  one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                equilibrium_populations(i,j,2) =  one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                equilibrium_populations(i,j,3) =  one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                equilibrium_populations(i,j,4) =  one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                equilibrium_populations(i,j,5) = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                equilibrium_populations(i,j,6) = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                equilibrium_populations(i,j,7) = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                equilibrium_populations(i,j,8) = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
            }            
        }


        /**
         * @brief Static method to update populations in D2Q9
         * 
         * @param populations tensor containing the populations of the nodes
         * @param fluid_nodes vector containing all the fluid nodes
         * @param global_rho tensor containing rho of all nodes
         * @param global_u tensor containing u of all nodes
         */
        static void update_macro(const Tensor<double, 3> &populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u)
        {
            Eigen::Index j, i;
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double rho, rhoinv;
            double ux, uy;

            for(size_t fnode = 0; fnode < fluid_nodes.size(); fnode++)
            {
                i = fluid_nodes[fnode].coords[0];
                j = fluid_nodes[fnode].coords[1];

                p0 = populations(i,j,0);
                p1 = populations(i,j,1);
                p2 = populations(i,j,2);
                p3 = populations(i,j,3);
                p4 = populations(i,j,4);
                p5 = populations(i,j,5);
                p6 = populations(i,j,6);
                p7 = populations(i,j,7);
                p8 = populations(i,j,8);

                rho = p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8;
                rhoinv = 1.0/rho;

                ux = rhoinv * (p1 + p5 + p8 - (p3 + p6 + p7));
                uy = rhoinv * (p2 + p5 + p6 - (p4 + p7 + p8));

                global_rho(i,j) = rho;  
                global_u(i,j,0) = ux;
                global_u(i,j,1) = uy;
            }
        }
    };

    template<std::size_t dim>
    class OMPDefaultEquilibrium : OMPTag
    {}; 

    template<>
    class OMPDefaultEquilibrium<2> : OMPTag {    
    
    public:
    static void calc_equilibrium(std::vector<Point<2>> &fluid_nodes, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &global_u, Tensor<double,2> &global_rho){
        double rho,ux,uy,u2;
        constexpr double two9 = (2.0/9.0);
        constexpr double one18 = (1.0/18.0);
        constexpr double one36 = (1.0/36.0);
        
        #pragma omp parallel for private(rho,ux,uy,u2)
        for (size_t k = 0; k < fluid_nodes.size(); k++)
        {
            size_t i = fluid_nodes[k].coords[0];
            size_t j = fluid_nodes[k].coords[1];
            rho = global_rho(i,j);
            ux = global_u(i,j,0);
            uy = global_u(i,j,1);
            u2 = ux * ux + uy * uy;
            equilibrium_populations(i,j,0) =  two9 * rho * (2 - 3*u2);
            equilibrium_populations(i,j,1) =  one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
            equilibrium_populations(i,j,2) =  one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
            equilibrium_populations(i,j,3) =  one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
            equilibrium_populations(i,j,4) =  one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
            equilibrium_populations(i,j,5) = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
            equilibrium_populations(i,j,6) = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
            equilibrium_populations(i,j,7) = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
            equilibrium_populations(i,j,8) = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
        }            
    }

    static void calc_equilibrium(std::vector<BoundaryPoint<2>> &fluid_nodes, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &global_u, Tensor<double,2> &global_rho){
        double rho,ux,uy,u2;
        constexpr double two9 = (2.0/9.0);
        constexpr double one18 = (1.0/18.0);
        constexpr double one36 = (1.0/36.0);
        
        #pragma omp parallel for private(rho,ux,uy,u2)
        for (size_t k = 0; k < fluid_nodes.size(); k++)
        {
            size_t i = fluid_nodes[k].coords[0];
            size_t j = fluid_nodes[k].coords[1];
            rho = global_rho(i,j);
            ux = global_u(i,j,0);
            uy = global_u(i,j,1);
            u2 = ux * ux + uy * uy;
            equilibrium_populations(i,j,0) =  two9 * rho * (2 - 3*u2);
            equilibrium_populations(i,j,1) =  one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
            equilibrium_populations(i,j,2) =  one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
            equilibrium_populations(i,j,3) =  one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
            equilibrium_populations(i,j,4) =  one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
            equilibrium_populations(i,j,5) = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
            equilibrium_populations(i,j,6) = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
            equilibrium_populations(i,j,7) = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
            equilibrium_populations(i,j,8) = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
        }            
    }


    /**
     * @brief Static method to update populations in D2Q9
     * 
     * @param populations tensor containing the populations of the nodes
     * @param fluid_nodes vector containing all the fluid nodes
     * @param global_rho tensor containing rho of all nodes
     * @param global_u tensor containing u of all nodes
     */
    static void update_macro(const Tensor<double, 3> &populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u)
    {
        Eigen::Index j, i;
        double p0, p1, p2, p3, p4, p5, p6, p7, p8;
        double rho, rhoinv;
        double ux, uy;

        #pragma omp parallel for private(i, j, p0, p1, p2, p3, p4, p5, p6, p7, p8, rho, rhoinv, ux, uy)
        for(size_t fnode = 0; fnode < fluid_nodes.size(); fnode++)
        {
            i = fluid_nodes[fnode].coords[0];
            j = fluid_nodes[fnode].coords[1];

            p0 = populations(i,j,0);
            p1 = populations(i,j,1);
            p2 = populations(i,j,2);
            p3 = populations(i,j,3);
            p4 = populations(i,j,4);
            p5 = populations(i,j,5);
            p6 = populations(i,j,6);
            p7 = populations(i,j,7);
            p8 = populations(i,j,8);

            rho = p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8;
            rhoinv = 1.0/rho;

            ux = rhoinv * (p1 + p5 + p8 - (p3 + p6 + p7));
            uy = rhoinv * (p2 + p5 + p6 - (p4 + p7 + p8));

            global_rho(i,j) = rho;  
            global_u(i,j,0) = ux;
            global_u(i,j,1) = uy;
        }
    }

};

    template<
        std::size_t dim
        >
    class STDExecDefaultEquilibrium : STDExecTag
    {};

    template<>
    class STDExecDefaultEquilibrium<2> : STDExecTag {    
    
    public:
        static void calc_equilibrium(std::vector<Point<2>> &fluid_nodes, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &global_u, Tensor<double,2> &global_rho){
            double rho,ux,uy,u2;
            constexpr double two9 = (2.0/9.0);
            constexpr double one18 = (1.0/18.0);
            constexpr double one36 = (1.0/36.0);
            
            //for (size_t k = 0; k < fluid_nodes.size(); k++)
            std::for_each(std::execution::par,fluid_nodes.begin(),fluid_nodes.end(),[&](const auto fluid_node)
            {
                size_t i = fluid_node.coords[0];
                size_t j = fluid_node.coords[1];
                rho = global_rho(i,j);
                ux = global_u(i,j,0);
                uy = global_u(i,j,1);
                u2 = ux * ux + uy * uy;
                equilibrium_populations(i,j,0) =  two9 * rho * (2 - 3*u2);
                equilibrium_populations(i,j,1) =  one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                equilibrium_populations(i,j,2) =  one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                equilibrium_populations(i,j,3) =  one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                equilibrium_populations(i,j,4) =  one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                equilibrium_populations(i,j,5) = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                equilibrium_populations(i,j,6) = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                equilibrium_populations(i,j,7) = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                equilibrium_populations(i,j,8) = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
            });            
        }

        static void calc_equilibrium(std::vector<BoundaryPoint<2>> &fluid_nodes, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &global_u, Tensor<double,2> &global_rho){
            double rho,ux,uy,u2;
            constexpr double two9 = (2.0/9.0);
            constexpr double one18 = (1.0/18.0);
            constexpr double one36 = (1.0/36.0);
            
            std::for_each(std::execution::par,fluid_nodes.begin(),fluid_nodes.end(),[&](const auto fluid_node)
            {
                size_t i = fluid_node.coords[0];
                size_t j = fluid_node.coords[1];
                rho = global_rho(i,j);
                ux = global_u(i,j,0);
                uy = global_u(i,j,1);
                u2 = ux * ux + uy * uy;
                equilibrium_populations(i,j,0) =  two9 * rho * (2 - 3*u2);
                equilibrium_populations(i,j,1) =  one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                equilibrium_populations(i,j,2) =  one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                equilibrium_populations(i,j,3) =  one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                equilibrium_populations(i,j,4) =  one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                equilibrium_populations(i,j,5) = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                equilibrium_populations(i,j,6) = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                equilibrium_populations(i,j,7) = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                equilibrium_populations(i,j,8) = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
            });            
        }


        /**
         * @brief Static method to update populations in D2Q9
         * 
         * @param populations tensor containing the populations of the nodes
         * @param fluid_nodes vector containing all the fluid nodes
         * @param global_rho tensor containing rho of all nodes
         * @param global_u tensor containing u of all nodes
         */
        static void update_macro(const Tensor<double, 3> &populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u)
        {
            Eigen::Index j, i;
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double rho, rhoinv;
            double ux, uy;

            //for(size_t fnode = 0; fnode < fluid_nodes.size(); fnode++)
            std::for_each(std::execution::par,fluid_nodes.begin(),fluid_nodes.end(),[&](const auto fluid_node)
            {
                i = fluid_node.coords[0];
                j = fluid_node.coords[1];

                p0 = populations(i,j,0);
                p1 = populations(i,j,1);
                p2 = populations(i,j,2);
                p3 = populations(i,j,3);
                p4 = populations(i,j,4);
                p5 = populations(i,j,5);
                p6 = populations(i,j,6);
                p7 = populations(i,j,7);
                p8 = populations(i,j,8);

                rho = p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8;
                rhoinv = 1.0/rho;

                ux = rhoinv * (p1 + p5 + p8 - (p3 + p6 + p7));
                uy = rhoinv * (p2 + p5 + p6 - (p4 + p7 + p8));

                global_rho(i,j) = rho;  
                global_u(i,j,0) = ux;
                global_u(i,j,1) = uy;
            });
        }
    };

    template<
        std::size_t dim
        >
    class OpenACCDefaultEquilibrium : OpenACCTag
    {};

    template<>
    class OpenACCDefaultEquilibrium<2> : OpenACCTag {    
    
    public:
        static void calc_equilibrium(std::vector<Point<2>> &fluid_nodes, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &global_u, Tensor<double,2> &global_rho){
            double rho,ux,uy,u2;
            constexpr double two9 = (2.0/9.0);
            constexpr double one18 = (1.0/18.0);
            constexpr double one36 = (1.0/36.0);

            auto n_rows = equilibrium_populations.dimensions()[0];
            auto n_cols = equilibrium_populations.dimensions()[1];
            const std::size_t num_directions = 9;

            double* equilibrium_populations_buffer = new double[n_rows * n_cols * num_directions];

            for (Eigen::Index i = 0; i < n_rows; ++i) {
                for (Eigen::Index j = 0; j < n_cols; ++j) {
                    for (std::size_t d = 0; d < num_directions; ++d) {
                        equilibrium_populations_buffer[i * n_cols * num_directions + j * num_directions + d] = equilibrium_populations(i, j, d);
                    }
                }
            }

            #pragma acc data copy(equilibrium_populations_buffer[0:n_rows * n_cols * num_directions])
            #pragma acc cache(equilibrium_populations_buffer[0:n_rows * n_cols * num_directions])
            {
                #pragma acc parallel loop
                for (std::size_t k = 0; k < fluid_nodes.size(); k++)
                {
                    Eigen::Index i = fluid_nodes[k].coords[0];
                    Eigen::Index j = fluid_nodes[k].coords[1];
                    double rho = global_rho(i, j);
                    double ux = global_u(i, j, 0);
                    double uy = global_u(i, j, 1);
                    double u2 = ux * ux + uy * uy;

                    std::size_t base_idx = i * n_cols * num_directions + j * num_directions;

                    equilibrium_populations_buffer[base_idx + 0] =  two9 * rho * (2 - 3 * u2);
                    equilibrium_populations_buffer[base_idx + 1] =  one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                    equilibrium_populations_buffer[base_idx + 2] =  one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                    equilibrium_populations_buffer[base_idx + 3] =  one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                    equilibrium_populations_buffer[base_idx + 4] =  one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                    equilibrium_populations_buffer[base_idx + 5] = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    equilibrium_populations_buffer[base_idx + 6] = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    equilibrium_populations_buffer[base_idx + 7] = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    equilibrium_populations_buffer[base_idx + 8] = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                }
            }

            for (Eigen::Index i = 0; i < n_rows; ++i) {
                for (Eigen::Index j = 0; j < n_cols; ++j) {
                    for (std::size_t d = 0; d < num_directions; ++d) {
                        equilibrium_populations(i, j, d) = equilibrium_populations_buffer[i * n_cols * num_directions + j * num_directions + d];
                    }
                }
            }

            delete[] equilibrium_populations_buffer;
        }

        static void calc_equilibrium(std::vector<BoundaryPoint<2>> &fluid_nodes, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &global_u, Tensor<double,2> &global_rho){
            double rho,ux,uy,u2;
            constexpr double two9 = (2.0/9.0);
            constexpr double one18 = (1.0/18.0);
            constexpr double one36 = (1.0/36.0);
            
            auto n_rows = equilibrium_populations.dimensions()[0];
            auto n_cols = equilibrium_populations.dimensions()[1];
            const std::size_t num_directions = 9;

            double* equilibrium_populations_buffer = new double[n_rows * n_cols * num_directions];

            for (Eigen::Index i = 0; i < n_rows; ++i) {
                for (Eigen::Index j = 0; j < n_cols; ++j) {
                    for (std::size_t d = 0; d < num_directions; ++d) {
                        equilibrium_populations_buffer[i * n_cols * num_directions + j * num_directions + d] = equilibrium_populations(i, j, d);
                    }
                }
            }

            #pragma acc data copy(equilibrium_populations_buffer[0:n_rows * n_cols * num_directions])
            #pragma acc cache(equilibrium_populations_buffer[0:n_rows * n_cols * num_directions])
            {
                #pragma acc parallel loop
                for (std::size_t k = 0; k < fluid_nodes.size(); k++)
                {
                    Eigen::Index i = fluid_nodes[k].coords[0];
                    Eigen::Index j = fluid_nodes[k].coords[1];
                    double rho = global_rho(i, j);
                    double ux = global_u(i, j, 0);
                    double uy = global_u(i, j, 1);
                    double u2 = ux * ux + uy * uy;

                    std::size_t base_idx = i * n_cols * num_directions + j * num_directions;

                    equilibrium_populations_buffer[base_idx + 0] =  two9 * rho * (2 - 3 * u2);
                    equilibrium_populations_buffer[base_idx + 1] =  one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                    equilibrium_populations_buffer[base_idx + 2] =  one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                    equilibrium_populations_buffer[base_idx + 3] =  one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                    equilibrium_populations_buffer[base_idx + 4] =  one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                    equilibrium_populations_buffer[base_idx + 5] = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    equilibrium_populations_buffer[base_idx + 6] = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    equilibrium_populations_buffer[base_idx + 7] = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    equilibrium_populations_buffer[base_idx + 8] = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                }
            }

            for (Eigen::Index i = 0; i < n_rows; ++i) {
                for (Eigen::Index j = 0; j < n_cols; ++j) {
                    for (std::size_t d = 0; d < num_directions; ++d) {
                        equilibrium_populations(i, j, d) = equilibrium_populations_buffer[i * n_cols * num_directions + j * num_directions + d];
                    }
                }
            }

            delete[] equilibrium_populations_buffer;

        }


        /**
         * @brief Static method to update populations in D2Q9
         * 
         * @param populations tensor containing the populations of the nodes
         * @param fluid_nodes vector containing all the fluid nodes
         * @param global_rho tensor containing rho of all nodes
         * @param global_u tensor containing u of all nodes
         */
        static void update_macro(const Tensor<double, 3> &populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u)
        {
            Eigen::Index j, i;
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double rho, rhoinv;
            double ux, uy;

            #pragma acc parallel loop
            for(size_t fnode = 0; fnode < fluid_nodes.size(); fnode++)
            {
                i = fluid_nodes[fnode].coords[0];
                j = fluid_nodes[fnode].coords[1];

                p0 = populations(i,j,0);
                p1 = populations(i,j,1);
                p2 = populations(i,j,2);
                p3 = populations(i,j,3);
                p4 = populations(i,j,4);
                p5 = populations(i,j,5);
                p6 = populations(i,j,6);
                p7 = populations(i,j,7);
                p8 = populations(i,j,8);

                rho = p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8;
                rhoinv = 1.0/rho;

                ux = rhoinv * (p1 + p5 + p8 - (p3 + p6 + p7));
                uy = rhoinv * (p2 + p5 + p6 - (p4 + p7 + p8));

                global_rho(i,j) = rho;  
                global_u(i,j,0) = ux;
                global_u(i,j,1) = uy;
            }
        }
    };



} // namespace llalbm::core::equilibrium

#endif