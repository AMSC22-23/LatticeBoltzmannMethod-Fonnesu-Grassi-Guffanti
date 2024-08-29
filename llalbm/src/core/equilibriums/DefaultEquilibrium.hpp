#ifndef LLALBM_DEFAULT_EQUILIBRIUM_HPP
#define LLALBM_DEFAULT_EQUILIBRIUM_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
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
} // namespace llalbm::core::equilibrium

#endif