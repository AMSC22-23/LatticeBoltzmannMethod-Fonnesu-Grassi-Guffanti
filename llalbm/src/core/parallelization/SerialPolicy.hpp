/**
 * @file SerialPolicy.hpp
 * @author Luca Guffanti
 * @brief Contains the Serial parallelization policy, which is the base one 
 * 
 */

#ifndef LLALBM_SERIALPOLICY_HPP
#define LLALBM_SERIALPOLICY_HPP


// =========== LLALBM INCLUDES ========
#include "../PolicyTypes.hpp"
#include "../../utils/aliases.hpp"
// ====================================

// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

namespace llalbm::core
{
    /**
     * @brief Serial policy, which is the base one.
     * 
     */
    template <std::size_t dim, typename Configuration>
    requires IsGloballySerial<Configuration>
    class SerialPolicy : public SequentialTag {
    public:
        SerialPolicy() {
            std::cout << "Generic Serial Policy is not implemented" << std::endl;
        }
    };


    template < typename Configuration >
    class SerialPolicy<2, Configuration> : public SequentialTag {
    public:

        using CollisionPolicy = typename Configuration::collision_policy_t;
        using WallPolicy = typename Configuration::wall_policy_t;
        using ObstaclePolicy = typename Configuration::obstacle_policy_t;
        using InletPolicy = typename Configuration::inlet_policy_t;
        using OutletPolicy = typename Configuration::outlet_policy_t;
        using InitializationPolicy = typename Configuration::initialization_policy_t;
        using EquilibriumPolicy = typename Configuration::equilibrium_policy_t;


        // ========================================================================================= 
        //  UPDATE OF THE BOUNDARIES 
        // ========================================================================================= 
        static void update_inlet_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
        {
            InletPolicy::update_boundaries(populations, boundary_coord, global_rho, global_u);
            return;
        }

        static void update_outlet_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
        {
            OutletPolicy::update_boundaries(populations, boundary_coord, global_rho, global_u);
            return;
        }

        static void update_domain_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
        {
            WallPolicy::update_boundaries(populations, boundary_coord, global_rho, global_u);
            return;
        }

        static void update_obstacle_boundaries(Tensor<double, 3> &populations, std::vector<ObstaclePoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
        {
            ObstaclePolicy::update_boundaries(populations, boundary_coord, global_rho, global_u);
            return;
        }

        // ========================================================================================= 
        // ========================================================================================= 
        // ========================================================================================= 

        // ========================================================================================= 
        //  COLLISION AND STREAMING
        // ========================================================================================= 

        static void collide(const Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u, const double time_step)
        {
            CollisionPolicy::collide(populations, equilibrium_populations, after_collision_populations, fluid_nodes, global_rho, global_u, time_step);
            return;
        }
        static void stream(Tensor<double, 3> &populations, Tensor<double, 3> &after_collision_populations, const std::vector<Point<2>> &fluid_nodes)
        {   
            CollisionPolicy::stream(populations, after_collision_populations, fluid_nodes);
            return;
        }   

        static void collide_open_inlet_boundary(const Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, const std::vector<BoundaryPoint<2>> &open_boundary_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u, const double time_step)
        {
            CollisionPolicy::collide_open_boundary(populations, equilibrium_populations, after_collision_populations, open_boundary_nodes, global_rho, global_u, time_step);
            return;
        }

        static void collide_open_outlet_boundary(const Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, const std::vector<BoundaryPoint<2>> &open_boundary_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u, const double time_step)
        {
            CollisionPolicy::collide_open_boundary(populations, equilibrium_populations, after_collision_populations, open_boundary_nodes, global_rho, global_u, time_step);
            return;
        }

        static void stream_open_inlet_boundary(Tensor<double, 3> &populations, Tensor<double, 3> &after_collision_populations, const std::vector<BoundaryPoint<2>> &open_boundary_nodes)
        {
            CollisionPolicy::stream_open_boundary(populations, after_collision_populations, open_boundary_nodes);
            return;
        }

        static void stream_open_outlet_boundary(Tensor<double, 3> &populations, Tensor<double, 3> &after_collision_populations, const std::vector<BoundaryPoint<2>> &open_boundary_nodes)
        {
            CollisionPolicy::stream_open_boundary(populations, after_collision_populations, open_boundary_nodes);
            return;
        }

        // ========================================================================================= 
        // ========================================================================================= 
        // ========================================================================================= 


        // ========================================================================================= 
        //  INITIALIZATION OF THE DOMAIN AND EXTERNAL UPDATES 
        // ========================================================================================= 


        static void attach_nodes( std::vector<BoundaryPoint<2>>& inlet_nodes_, std::vector<BoundaryPoint<2>>& outlet_nodes_)
        {
            InitializationPolicy::attach_nodes(inlet_nodes_, outlet_nodes_);
            return;
        }

        static void update_macro(const Tensor<double, 3> &populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u)
        {
            InitializationPolicy::update_macro(populations, fluid_nodes, global_rho, global_u);
            return;
        }

        static void update_nodes(const double& time_step, Eigen::Tensor<double, 3>& velocity_tensor, Eigen::Tensor<double, 2> density_tensor)
        {
            InitializationPolicy::update_nodes(time_step, velocity_tensor, density_tensor);
            return;
        }
        // ========================================================================================= 
        // ========================================================================================= 
        // ========================================================================================= 


        // ========================================================================================= 
        //  COMPUTATION OF EQUILIBRIUM 
        // ========================================================================================= 

        static void calc_equilibrium(std::vector<Point<2>> &fluid_nodes, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &global_u, Tensor<double,2> &global_rho)
        {
            EquilibriumPolicy::calc_equilibrium(fluid_nodes, equilibrium_populations, global_u, global_rho);
            return;
        }

        static void calc_equilibrium(std::vector<BoundaryPoint<2>> &fluid_nodes, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &global_u, Tensor<double,2> &global_rho)
        {
            EquilibriumPolicy::calc_equilibrium(fluid_nodes, equilibrium_populations, global_u, global_rho);
            return;
        }
    };
} // namespace llalbm::core


#endif // LLALBM_SERIALPOLICY_HPP