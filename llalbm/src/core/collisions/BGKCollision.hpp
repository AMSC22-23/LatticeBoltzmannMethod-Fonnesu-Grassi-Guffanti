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
#include "Logger.hpp"
#include "aliases.hpp" 
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
    public:
        static double tau;
        const double tauinv = 1.0/tau;
        const double omtauinv = 1.0 - tauinv;
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
        }

        /**
         * @brief Calculates collisions using bgk
         * 
         * @param populations of the nodes
         * @param equilibrium_populations of the nodes
         * @param after_collision_populations to be streamed
         */
        void stream_collide(Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, std::vector<Point<2>> fluid_nodes, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
        {
            for(size_t fnode = 0; fnode < fluid_nodes.size(); fnode++)
            {
                double x = fluid_nodes[fnode][0];
                double y = fluid_nodes[fnode][1];

                double p0 = populations(x, y, 0);
                double p1 = populations(x-1, y, 1);
                double p2 = populations(x, y-1, 2);
                double p3 = populations(x+1, y, 3);
                double p4 = populations(x, y+1, 4);
                double p5 = populations(x-1, y-1, 5);
                double p6 = populations(x+1, y-1, 6);
                double p7 = populations(x+1, y+1, 7);
                double p8 = populations(x-1, y+1, 8);

                double rho = p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8;
                double rhoinv = 1.0/rho;

                double ux = rhoinv * (p1 + p5 + p8 -(p3 + p6 + p7));
                double uy = rhoinv * (p2 + p5 + p6 -(p4+p7+p8));

                double tw0r = tauinv * D2Q9(0, 2) * rho;
                double twsr = tauinv * D2Q9(1, 2) * rho;
                double twdr = tauinv * D2Q9(5, 2) * rho;
            
                double omusq = 1.0 - 1.5 * (ux * ux + uy*uy);

                double tux = 3.0*ux;
                double tuy = 3.0*uy;
                after_collision_populations(x,y,0)= omtauinv*p0 + tw0r*(omusq);

                double cidot3u = tux; //(?)
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
        }
    };
    // initialization of the relaxation constant in the 2-D BGK collision operator.
    //double BGKCollisionPolicy<2>::tau = 0.0;

}; // namespace llalbm::core::collisions
