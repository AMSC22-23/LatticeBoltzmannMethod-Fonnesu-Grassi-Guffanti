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
        void collide(Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, std::vector<Point<2>> fluid_nodes)
        {

        }
    };
    // initialization of the relaxation constant in the 2-D BGK collision operator.
    double BGKCollisionPolicy<2>::tau = 0.0;

}; // namespace llalbm::core::collisions