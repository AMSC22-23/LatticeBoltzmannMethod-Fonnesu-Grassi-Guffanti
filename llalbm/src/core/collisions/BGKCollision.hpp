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
        Eigen::Index x, y;
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
        void stream_collide(const Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u, bool save)
        {

            for(size_t fnode = 0; fnode < fluid_nodes.size(); fnode++)
            {
                x = fluid_nodes[fnode].coords[0];
                y = fluid_nodes[fnode].coords[1];

                p0 = populations(x, y, 0);
                p1 = populations(x-1, y, 1);
                p2 = populations(x, y-1, 2);
                p3 = populations(x+1, y, 3);
                p4 = populations(x, y+1, 4);
                p5 = populations(x-1, y-1, 5);
                p6 = populations(x+1, y-1, 6);
                p7 = populations(x+1, y+1, 7);
                p8 = populations(x-1, y+1, 8);

                rho = p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8;
                rhoinv = 1.0/rho;

                ux = rhoinv * (p1 + p5 + p8 -(p3 + p6 + p7));
                uy = rhoinv * (p2 + p5 + p6 -(p4+p7+p8));

                if(save)
                {
                    global_rho(x,y) = rho;
                    global_u(x,y,0) = ux;
                    global_u(x,y,1) = uy;
                }

                tw0r = tauinv * D2Q9(0, 2) * rho;
                twsr = tauinv * D2Q9(1, 2) * rho;
                twdr = tauinv * D2Q9(5, 2) * rho;
            
                omusq = 1.0 - 1.5 * (ux * ux + uy*uy);

                tux = 3.0*ux;
                tuy = 3.0*uy;
                
                after_collision_populations(x,y,0)= omtauinv*p0 + tw0r*(omusq);

                cidot3u = tux; //(?)
                after_collision_populations(x,y,1)= omtauinv*p1 + twsr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
                after_collision_populations(x,y,2)= omtauinv*p2 + twsr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
                after_collision_populations(x,y,3)= omtauinv*p3 + twsr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
                after_collision_populations(x,y,4)= omtauinv*p4 + twsr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));

                cidot3u = tux + tuy;
                after_collision_populations(x,y,5)= omtauinv*p5 + twdr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
                after_collision_populations(x,y,6)= omtauinv*p6 + twdr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
                after_collision_populations(x,y,7)= omtauinv*p7 + twdr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
                after_collision_populations(x,y,8)= omtauinv*p8 + twdr *(omusq + cidot3u * (1.0 + 0.5*cidot3u));
            }

            // TODO: streaming
            /* Dopo che sono state calcolate le after_collision_populations di TUTTI i fluid nodes, bisogna fare lo streaming.
            3 modi (aggiungere altri in caso):
            - se possibile aggiungere dentro il loop di stream_collide (post-sincronizzazione) per evitare di averne un altro;
            - aggiungere un nuovo loop in stream_collide;
            - creare un nuovo metodo per lo streaming.
             */
        }
    };
    // initialization of the relaxation constant in the 2-D BGK collision operator.
    double BGKCollisionPolicy<2>::tau = 0.0;
    double BGKCollisionPolicy<2>::tauinv = 0.0;
    double BGKCollisionPolicy<2>::omtauinv = 0.0;

}; // namespace llalbm::core::collisions
