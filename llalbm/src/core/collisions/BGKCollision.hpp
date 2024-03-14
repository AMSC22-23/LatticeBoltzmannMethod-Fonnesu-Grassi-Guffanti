/**
 * @file BGKCollision.hpp File containing the templates for the BGK collision policy
 * @brief The BGK collision in the Lattice-Boltzmann method is among the 
 * 
 */

// =========== STL INCLUDES ===========
#include <vector>
#include <map>
#include <array>
#include <cassert>
// ======================================

// =========== EIGEN INCLUDES ===========
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "../../utils/loggers/Logger.hpp"
#include "../../utils/aliases.hpp" 
// =======================================

namespace llalbm::core::collisions 
{

    using namespace llalbm::util::logger;

    /**
     * @brief Most generic version of the BGK collision operator.
     * 
     * @tparam dim the total number of spatial dimensions.
     */
    template<
        std::size_t dim    
    >
    class BGKCollisionPolicy
    {
    private:
        /// @brief Logger of the collision policy
        Logger l;

    public:
        BGKCollisionPolicy() 
        : l ("GENERIC BGK", std::cout)
        {
            l.error("GENERIC BGK COLLISION NOT IMPLEMENTED");
            assert(false);
        }
    };


    /**
     * @brief Specialization of the BGK collision operator for the 
     * use-case of two spatial dimensions.
     */
    template<>
    class BGKCollisionPolicy<2>
    {
    private:
        //std::array<Eigen::Index, 2> lattice_dimensions;
        Eigen::Index j, i;
        double p0, p1, p2, p3, p4, p5, p6, p7, p8;
        double rho, rhoinv;
        double ux, uy;
        double tw0r, twsr, twdr;
        double omusq, tux, tuy, cidot3u;
    public:
        static double tau;
        static double tauinv;
        static double omtauinv;
        //std::array<Eigen::Index, 2> lattice_nodes; 
        
        /**
         * @brief Initializes the BGK collision policy with the
         * relaxation constant.
         * 
         * @param tau_ the relaxation constant
         */
        static void initialize(double tau_)
        {
            tau = tau_;
            tauinv = 1.0/tau;
            omtauinv = 1.0 - tauinv;
        }

        /**
         * @brief method that calculates the new populations of a node, its rho and its velocities
         * 
         * @param populations of the node
         * @param equilibrium_populations of the node
         * @param after_collision_populations of the node to be streamed
         * @param fluid_nodes vector containing coordinates of all the fluid nodes of the lattice
         * @param global_rho tensor containing all densities
         * @param global_u tensor containing all velocities 
         * @param save boolean used to save rho and u when necessary
         */
        void collide(const Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u)
        {

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

                tw0r = tauinv * D2Q9(0, 2) * rho;
                twsr = tauinv * D2Q9(1, 2) * rho;
                twdr = tauinv * D2Q9(5, 2) * rho;
            
                omusq = 1.0 - 1.5 * (ux * ux + uy*uy);

                tux = 3.0*ux;
                tuy = 3.0*uy;
                
                after_collision_populations(i,j,0)= omtauinv*p0 + tw0r*(omusq);

                cidot3u = tux;
                after_collision_populations(i,j,1)= omtauinv*p1 + twsr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
                after_collision_populations(i,j,2)= omtauinv*p2 + twsr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
                after_collision_populations(i,j,3)= omtauinv*p3 + twsr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
                after_collision_populations(i,j,4)= omtauinv*p4 + twsr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));

                cidot3u = tux + tuy;
                after_collision_populations(i,j,5)= omtauinv*p5 + twdr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
                after_collision_populations(i,j,6)= omtauinv*p6 + twdr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
                after_collision_populations(i,j,7)= omtauinv*p7 + twdr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
                after_collision_populations(i,j,8)= omtauinv*p8 + twdr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
            }
        }

        /**
         * @brief method that performs the streaming from after_collision_populations to populations for every fluid node
         * 
         * @param populations of fluid nodes to be updated
         * @param after_collision_populations of fluid nodes to be streamed
         * @param fluid_nodes vector containing coordinates of all the fluid nodes of the lattice
         */
        void stream(Tensor<double, 3> &populations, Tensor<double, 3> &after_collision_populations, const std::vector<Point<2>> &fluid_nodes)
        {
            for(size_t fnode = 0; fnode < fluid_nodes.size(); fnode++)
            {
                i = fluid_nodes[fnode].coords[0];
                j = fluid_nodes[fnode].coords[1];

                populations(i, j, 0) = after_collision_populations(i,j,0);
                populations(i, j-1, 1) = after_collision_populations(i,j,1);
                populations(i+1, j, 2) = after_collision_populations(i,j,2);
                populations(i, j+1, 3) = after_collision_populations(i,j,3);
                populations(i-1, j, 4) = after_collision_populations(i,j,4);
                populations(i+1, j-1, 5) = after_collision_populations(i,j,5);
                populations(i+1, j+1, 6) = after_collision_populations(i,j,6);
                populations(i-1, j+1, 7) = after_collision_populations(i,j,7);
                populations(i-1, j-1, 8) = after_collision_populations(i,j,8);
            }
        }
    };

    // initialization of the relaxation constant in the 2-D BGK collision operator.
    double BGKCollisionPolicy<2>::tau = 0.0;
    double BGKCollisionPolicy<2>::tauinv = 0.0;
    double BGKCollisionPolicy<2>::omtauinv = 0.0;

}; // namespace llalbm::core::collisions
