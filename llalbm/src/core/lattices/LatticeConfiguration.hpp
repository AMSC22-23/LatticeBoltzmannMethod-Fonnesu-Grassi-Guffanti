 /**
 * @file LatticeConfiguration.hpp
 * @brief Contains the class that describes the configuration of a lattice,
 * in terms of dimensions and policies
 */

#ifndef LLALBM_LATTICECONFIGURATION_HPP
#define LLALBM_LATTICECONFIGURATION_HPP

// =========== STL INCLUDES ===========
#include <cstddef> // for std::size_t
// ====================================

// ========== LLALBM INCLUDES =========
#include "../PolicyTypes.hpp" // for std::size_t
#include "../parallelization/SerialPolicy.hpp" // for serial execution
// ====================================
namespace llalbm::core {

/**
 * @brief Lattice configuration static object used to group policies used during LBM simulations.
 * 
 * @tparam dim Physical dimensions of the lattice.
 * @tparam CollisionPolicy Policy managing the interaction between nodes.
 * @tparam WallPolicy Policy managing the interaction between a fluid node and a wall node.
 * @tparam ObstaclePolicy Policy managing the interaction between a fluid node and an internal obstacle.
 * @tparam InletPolicy Policy managing the evaluation of open boundary conditions for inlets.
 * @tparam OutletPolicy Policy managing the evaluation of open boundary conditions for outlets.
 * @tparam InitializationPolicy Policy managing the initialization of the inlets and outlets.
 * @tparam EquilibriumPolicy Policy managing the computation of the equilibrium distributions.
 */
template<
    std::size_t dim,                  // Spatial dimensions of the simulation.
    typename CollisionPolicy,         // Policy managing the interaction between two fluid nodes.
    typename WallPolicy,              // Policy managing the interaction between a fluid node and a wall node.
    typename ObstaclePolicy,          // Policy managing the interaction between a fluid node and an internal obstacle.
    typename InletPolicy,             // Policy managing the interaction between an internal fluid node and an inlet node.
    typename OutletPolicy,            // Policy managing the interaction between an internal fluid node and an outlet node.
    typename InitializationPolicy,    // Policy managing the initialization of the lattice.
    typename EquilibriumPolicy
> requires
    // Checking that the policy types are correct.
    IsCollisionPolicy<CollisionPolicy> &&
    IsBoundaryPolicy<WallPolicy> &&
    IsBoundaryPolicy<ObstaclePolicy> &&
    IsBoundaryPolicy<InletPolicy> &&
    IsBoundaryPolicy<OutletPolicy> &&
    IsInitializationPolicy<InitializationPolicy>

class LatticeConfiguration {
public:
    static constexpr std::size_t dimensions = dim; // Dimensions of the lattice.
    using collision_policy_t = CollisionPolicy;    // Type alias for the collision policy.
    using wall_policy_t = WallPolicy;              // Type alias for the wall policy.
    using obstacle_policy_t = ObstaclePolicy;      // Type alias for the obstacle policy.
    using inlet_policy_t = InletPolicy;            // Type alias for the inlet policy.
    using outlet_policy_t = OutletPolicy;          // Type alias for the outlet policy.
    using initialization_policy_t = InitializationPolicy; // Type alias for the initialization policy.
    using equilibrium_policy_t = EquilibriumPolicy; // Type alias for the equilibrium policy.
};

} // namespace llalbm::core

#endif // LLALBM_LATTICECONFIGURATION_HPP
