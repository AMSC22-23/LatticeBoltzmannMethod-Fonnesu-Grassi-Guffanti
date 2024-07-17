/**
 * @file PolicyTypes.hpp
 * @author Luca Guffanti
 * @brief This file defines classes used to tag the policies
 * to explicit their type and parallelized nature. The following types
 * are used in combination with the concept infrastructure to enforce compile
 * time type checking in template parameters, to ensure consistency when considering
 * both the parallelization of the policy and their actual type (Inlet, Outlet, etc.)
 * 
 * 
 * Templating and compile time type validity and coherence check are quite tricky but
 * vital in the context of the llalbm library, as it's header only.
 * The introduction of concepts in the last few ISO C++ standards has made it possible
 * by statically verifying that the types passed as template parameters are "marked" or
 * "tagged" in some way to express their nature. Policies are therefore tagged both in relation
 * to their type (i.e. Collision) and the way they are parallelized (i.e. parallelization with OpenMP).
 * 
 * Verification is done directly by the @link llalbm::core::LatticeConfiguration LatticeConfiguration @endlink class, which
 * encapsulates the policy choice for a specific lattice, using the concepts defined in this file.
 * 
 */

#ifndef LLALBM_POLICYTYPES_HPP
#define LLALBM_POLICYTYPES_HPP

// =========== STL INCLUDES ===========
#include <concepts>
#include <type_traits>  
// ======================================


namespace llalbm::core
{

    // =========== POLICY TYPES ===========
    class CollisionPolicyTag {}; 
    class BoundaryPolicyTag {};
    class InitializationPolicyTag {};

    // =========== PARALLELIZATION POLICIES ===========
    class OMPTag {};
    class SequentialTag {};
    class MPITag {};

    // =========== POLICY CONCEPTS ===========
    template<typename T>
    concept IsCollisionPolicy = std::is_base_of_v<CollisionPolicyTag, T>;

    template<typename T>
    concept IsBoundaryPolicy = std::is_base_of_v<BoundaryPolicyTag, T>;

    template<typename T>
    concept IsOMP = std::is_base_of_v<OMPTag, T>;

    template<typename T>
    concept IsSequential = std::is_base_of_v<SequentialTag, T>;

    template<typename T>
    concept IsMPI = std::is_base_of_v<MPITag, T>;

    template <typename T>
    concept IsInitializationPolicy = std::is_base_of_v<InitializationPolicyTag, T>;

    template<
        typename CollisionPolicy,
        typename WallPolicy,
        typename ObstaclePolicy,
        typename InletPolicy,
        typename OutletPolicy,
        typename InitializationPolicy,
        typename EquilibriumPolicy>
    concept IsGloballySerial = std::conjunction_v<
        std::is_base_of<SequentialTag, CollisionPolicy>,
        std::is_base_of<SequentialTag, WallPolicy>,
        std::is_base_of<SequentialTag, ObstaclePolicy>,
        std::is_base_of<SequentialTag, InletPolicy>,
        std::is_base_of<SequentialTag, OutletPolicy>,
        std::is_base_of<SequentialTag, InitializationPolicy>,
        std::is_base_of<SequentialTag, EquilibriumPolicy>>;

    template<
        typename CollisionPolicy,
        typename WallPolicy,
        typename ObstaclePolicy,
        typename InletPolicy,
        typename OutletPolicy,
        typename InitializationPolicy,
        typename EquilibriumPolicy>
    concept IsGloballyOpenMP = std::conjunction_v<
        std::is_base_of<OMPTag, CollisionPolicy>,
        std::is_base_of<OMPTag, WallPolicy>,
        std::is_base_of<OMPTag, ObstaclePolicy>,
        std::is_base_of<OMPTag, InletPolicy>,
        std::is_base_of<OMPTag, OutletPolicy>,
        std::is_base_of<OMPTag, InitializationPolicy>,
        std::is_base_of<SequentialTag, EquilibriumPolicy>>;

    template<
        typename CollisionPolicy,
        typename WallPolicy,
        typename ObstaclePolicy,
        typename InletPolicy,
        typename OutletPolicy,
        typename InitializationPolicy,
        typename EquilibriumPolicy>
    concept IsGloballyMPI = std::conjunction_v<
        std::is_base_of<MPITag, CollisionPolicy>,
        std::is_base_of<MPITag, WallPolicy>,
        std::is_base_of<MPITag, ObstaclePolicy>,
        std::is_base_of<MPITag, InletPolicy>,
        std::is_base_of<MPITag, OutletPolicy>,
        std::is_base_of<MPITag, InitializationPolicy>,
        std::is_base_of<SequentialTag, EquilibriumPolicy>>;

}; // llalbm::core

#endif // LLALBM_POLICYTYPES_HPP