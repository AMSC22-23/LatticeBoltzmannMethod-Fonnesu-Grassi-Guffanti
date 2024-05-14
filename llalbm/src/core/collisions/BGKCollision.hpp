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
        Eigen::Index j, i;
        double p0, p1, p2, p3, p4, p5, p6, p7, p8;
        double t;

    public:
        static double tau;
        
        /**
         * @brief Initializes the BGK collision policy with the
         * relaxation constant.
         * 
         * @param tau_ the relaxation constant
         */
        static void initialize(double tau_)
        {
            tau = tau_;
            std::cout << "Inizialization done" << std::endl;
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
        void collide(const Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u, const double time_step)
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

                t = time_step / tau;

                after_collision_populations(i,j,0) = p0 - t * (p0 - equilibrium_populations(i,j,0));
                after_collision_populations(i,j,1) = p1 - t * (p1 - equilibrium_populations(i,j,1));
                after_collision_populations(i,j,2) = p2 - t * (p2 - equilibrium_populations(i,j,2));
                after_collision_populations(i,j,3) = p3 - t * (p3 - equilibrium_populations(i,j,3));
                after_collision_populations(i,j,4) = p4 - t * (p4 - equilibrium_populations(i,j,4));
                after_collision_populations(i,j,5) = p5 - t * (p5 - equilibrium_populations(i,j,5));
                after_collision_populations(i,j,6) = p6 - t * (p6 - equilibrium_populations(i,j,6));
                after_collision_populations(i,j,7) = p7 - t * (p7 - equilibrium_populations(i,j,7));
                after_collision_populations(i,j,8) = p8 - t * (p8 - equilibrium_populations(i,j,8));
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

}; // namespace llalbm::core::collisions
