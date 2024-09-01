/**
 * @file TRTCollision.hpp File containing the templates for the TRT collision policy
 * @brief The TRT collision in the Lattice-Boltzmann method is a two-relaxation-time collision operator.
 * It is used to simulate the flow of a fluid with a higher degree of accuracy than the TRT collision operator, and 
 * requires two parameters, the relaxation constant for the shear viscosity and another free parameter.
 * 
 * 
 */

#ifndef LLALBM_TRT_COLLISION_HPP
#define LLALBM_TRT_COLLISION_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <map>
#include <array>
#include <cassert>
#include <execution>
// ======================================

// =========== EIGEN INCLUDES ===========
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "../../utils/loggers/Logger.hpp"
#include "../../utils/aliases.hpp" 
#include "../PolicyTypes.hpp"
// =======================================

namespace llalbm::core::collisions 
{

    using namespace llalbm::util::logger;

    /**
     * @brief Most generic version of the TRT collision operator.
     * 
     * @tparam dim the total number of spatial dimensions.
     */
    template<std::size_t dim>
    class TRTCollisionPolicy : public CollisionPolicyTag, public SequentialTag
    {
    private:
        /// @brief Logger of the collision policy
        Logger l;

    public:
        TRTCollisionPolicy() 
        : l ("GENERIC TRT", std::cout)
        {
            l.error("GENERIC TRT COLLISION NOT IMPLEMENTED");
            assert(false);
        }
    };


    /**
     * @brief Specialization of the TRT collision operator for the 
     * use-case of two spatial dimensions.
     */
    template<>
    class TRTCollisionPolicy<2> : public CollisionPolicyTag, public SequentialTag
    {

    public:
        static double tau_even;
        static double tau_odd;
        static double sound_velocity;

        /**
         * @brief Initializes the TRT collision policy with the
         * relaxation constants.
         * 
         * @param tau_even the relaxation constant associated to even moments
         * @param tau_odd the relaxation constant associated to odd moments
         * @param delta_t the time step of the simulation
         * 
         */
        static void initialize(const double& tau_even_, const double& tau_odd_, const double& sound_velocity_)
        {
            tau_even = tau_even_;
            tau_odd = tau_odd_;
            sound_velocity = sound_velocity_;

            std::cout << "Inizialization done:" << std::endl;
            std::cout << " Tau+ : " << tau_even << std::endl;
            std::cout << " Tau- : " << tau_odd << std::endl;
            std::cout << " Sound velocity: " << sound_velocity << std::endl;
        }

        /**
         * @brief Computes the kinematic shear viscosity of the fluid
         * using the parameters of the collision model.
         * 
         * @return the kinematic viscosity of the fluid
         */
        static double compute_kinematic_viscosity()
        {
            return (tau_even - 0.5) * sound_velocity * sound_velocity;
        }

        /**
         * @brief Computes the magic parameter of the model, which shows
         * how much the relaxation constants control stability and accuracy
         * 
         * @note Lambda = 1/12 cancels the third order spatial error
         * @note Lambda =  1/6  cancels the fourth order spatial error
         * @note Lambda =  1/4  provides stable simulations
         * 
         * @return the magic parameter of the model, Lambda 
         */
        static double compute_magic_parameter()
        {
            std::cout <<  "Magic Parameter: " << (tau_even - 0.5)*(tau_odd- 0.5) << std::endl;
            return (tau_even - 0.5)*(tau_odd- 0.5);
        }

        /**
         * @brief Sets the value of te free paramer for odd velocities by enforcing a magic parameer value,
         * which is used to control the stability and accuracy of the simulation, and is set to 1/4 by default.
         * 
         * @param magic_parameter magic parameter target
         */
        static void enforce_magic_parameter(const double& magic_parameter = 1.0/4.0)
        {
            tau_odd = (magic_parameter)*(1.0/(tau_even-0.5))+0.5;
            std::cout << "Magic Parameter forced to " << magic_parameter << std::endl;
            std::cout << " tau+ : " << tau_even << std::endl;
            std::cout << " tau- : " << tau_odd << std::endl;
        }

        /**
         * @brief Computes the collision of populations of fluid nodes
         * 
         * @param populations of the node
         * @param equilibrium_populations of the node
         * @param after_collision_populations of the node to be streamed
         * @param fluid_nodes vector containing coordinates of all the fluid nodes of the lattice
         * @param global_rho tensor containing all densities
         * @param global_u tensor containing all velocities 
         */
        static void collide(const Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u, const double time_step)
        {

            Eigen::Index j, i;

            for(const auto& fluid_node : fluid_nodes)
            {
                i = fluid_node.coords[0];
                j = fluid_node.coords[1];
                after_collision_populations(i,j,0) = populations(i,j,0) - 
                    (1.0/tau_even) * (populations(i,j,0) - equilibrium_populations(i,j,0));
                
                after_collision_populations(i,j,1) = populations(i,j,1) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,1) + populations(i,j,3) - equilibrium_populations(i,j,1) - equilibrium_populations(i,j,3)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,1) - populations(i,j,3) - equilibrium_populations(i,j,1) + equilibrium_populations(i,j,3));
                
                after_collision_populations(i,j,2) = populations(i,j,2) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,2) + populations(i,j,4) - equilibrium_populations(i,j,2) - equilibrium_populations(i,j,4)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,2) - populations(i,j,4) - equilibrium_populations(i,j,2) + equilibrium_populations(i,j,4));
                
                after_collision_populations(i,j,3) = populations(i,j,3) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,3) + populations(i,j,1) - equilibrium_populations(i,j,3) - equilibrium_populations(i,j,1)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,3) - populations(i,j,1) - equilibrium_populations(i,j,3) + equilibrium_populations(i,j,1));
                
                after_collision_populations(i,j,4) = populations(i,j,4) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,4) + populations(i,j,2) - equilibrium_populations(i,j,4) - equilibrium_populations(i,j,2)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,4) - populations(i,j,2) - equilibrium_populations(i,j,4) + equilibrium_populations(i,j,2));
                
                after_collision_populations(i,j,5) = populations(i,j,5) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,5) + populations(i,j,7) - equilibrium_populations(i,j,5) - equilibrium_populations(i,j,7)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,5) - populations(i,j,7) - equilibrium_populations(i,j,5) + equilibrium_populations(i,j,7));
                
                after_collision_populations(i,j,6) = populations(i,j,6) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,6) + populations(i,j,8) - equilibrium_populations(i,j,6) - equilibrium_populations(i,j,8)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,6) - populations(i,j,8) - equilibrium_populations(i,j,6) + equilibrium_populations(i,j,8));
                
                after_collision_populations(i,j,7) = populations(i,j,7) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,7) + populations(i,j,5) - equilibrium_populations(i,j,7) - equilibrium_populations(i,j,5)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,7) - populations(i,j,5) - equilibrium_populations(i,j,7) + equilibrium_populations(i,j,5));
                
                after_collision_populations(i,j,8) = populations(i,j,8) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,8) + populations(i,j,6) - equilibrium_populations(i,j,8) - equilibrium_populations(i,j,6)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,8) - populations(i,j,6) - equilibrium_populations(i,j,8) + equilibrium_populations(i,j,6));
                
            }
        }

        /**
         * @brief Performs the streaming from after_collision_populations to populations for every fluid node
         * 
         * @param populations of fluid nodes to be updated
         * @param after_collision_populations of fluid nodes to be streamed
         * @param fluid_nodes vector containing coordinates of all the fluid nodes of the lattice
         */
        static void stream(Tensor<double, 3> &populations, Tensor<double, 3> &after_collision_populations, const std::vector<Point<2>> &fluid_nodes)
        {

            Eigen::Index j, i;

            
            for(const auto& fluid_node : fluid_nodes)
            {
                i = fluid_node.coords[0];
                j = fluid_node.coords[1];
                
                populations(i, j, 0) = after_collision_populations(i,j,0);
                populations(i, j+1, 1) = after_collision_populations(i,j,1);
                populations(i-1, j, 2) = after_collision_populations(i,j,2);
                populations(i, j-1, 3) = after_collision_populations(i,j,3);
                populations(i+1, j, 4) = after_collision_populations(i,j,4);
                populations(i-1, j+1, 5) = after_collision_populations(i,j,5);
                populations(i-1, j-1, 6) = after_collision_populations(i,j,6);
                populations(i+1, j-1, 7) = after_collision_populations(i,j,7);
                populations(i+1, j+1, 8) = after_collision_populations(i,j,8);
            }
        }

    /**
     * @brief Computes the collision of populations at the open boundaries
     * 
     * @param populations of the node
     * @param equilibrium_populations of the node
     * @param after_collision_populations of the node to be streamed
     * @param open_boundary_nodes vector containing coordinates of all the border nodes of the lattice
     * @param global_rho tensor containing all densities
     * @param global_u tensor containing all velocities 
     */
    static void collide_open_boundary(const Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, const std::vector<BoundaryPoint<2>> &open_boundary_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u, const double time_step)
    {

        Eigen::Index j, i;

        for(const auto& open_boundary : open_boundary_nodes)
        {
            i = open_boundary.coords[0];
            j = open_boundary.coords[1];
            after_collision_populations(i,j,0) = populations(i,j,0) - 
                (1.0/tau_even) * (populations(i,j,0) - equilibrium_populations(i,j,0));
            
            after_collision_populations(i,j,1) = populations(i,j,1) - 
                0.5 * (1.0/tau_even) * (populations(i,j,1) + populations(i,j,3) - equilibrium_populations(i,j,1) - equilibrium_populations(i,j,3)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,1) - populations(i,j,3) - equilibrium_populations(i,j,1) + equilibrium_populations(i,j,3));
            
            after_collision_populations(i,j,2) = populations(i,j,2) - 
                0.5 * (1.0/tau_even) * (populations(i,j,2) + populations(i,j,4) - equilibrium_populations(i,j,2) - equilibrium_populations(i,j,4)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,2) - populations(i,j,4) - equilibrium_populations(i,j,2) + equilibrium_populations(i,j,4));
            
            after_collision_populations(i,j,3) = populations(i,j,3) - 
                0.5 * (1.0/tau_even) * (populations(i,j,3) + populations(i,j,1) - equilibrium_populations(i,j,3) - equilibrium_populations(i,j,1)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,3) - populations(i,j,1) - equilibrium_populations(i,j,3) + equilibrium_populations(i,j,1));
            
            after_collision_populations(i,j,4) = populations(i,j,4) - 
                0.5 * (1.0/tau_even) * (populations(i,j,4) + populations(i,j,2) - equilibrium_populations(i,j,4) - equilibrium_populations(i,j,2)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,4) - populations(i,j,2) - equilibrium_populations(i,j,4) + equilibrium_populations(i,j,2));
            
            after_collision_populations(i,j,5) = populations(i,j,5) - 
                0.5 * (1.0/tau_even) * (populations(i,j,5) + populations(i,j,7) - equilibrium_populations(i,j,5) - equilibrium_populations(i,j,7)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,5) - populations(i,j,7) - equilibrium_populations(i,j,5) + equilibrium_populations(i,j,7));
            
            after_collision_populations(i,j,6) = populations(i,j,6) - 
                0.5 * (1.0/tau_even) * (populations(i,j,6) + populations(i,j,8) - equilibrium_populations(i,j,6) - equilibrium_populations(i,j,8)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,6) - populations(i,j,8) - equilibrium_populations(i,j,6) + equilibrium_populations(i,j,8));
            
            after_collision_populations(i,j,7) = populations(i,j,7) - 
                0.5 * (1.0/tau_even) * (populations(i,j,7) + populations(i,j,5) - equilibrium_populations(i,j,7) - equilibrium_populations(i,j,5)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,7) - populations(i,j,5) - equilibrium_populations(i,j,7) + equilibrium_populations(i,j,5));
            
            after_collision_populations(i,j,8) = populations(i,j,8) - 
                0.5 * (1.0/tau_even) * (populations(i,j,8) + populations(i,j,6) - equilibrium_populations(i,j,8) - equilibrium_populations(i,j,6)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,8) - populations(i,j,6) - equilibrium_populations(i,j,8) + equilibrium_populations(i,j,6));
            
        }
    }

    /**
     * @brief Performs the streaming from after_collision_populations to populations for every open boundary node
     * 
     * @param populations of border nodes to be updated
     * @param after_collision_populations of border nodes to be streamed
     * @param open_boundary_nodes vector containing coordinates of all the border nodes of the lattice
     */
    static void stream_open_boundary(Tensor<double, 3> &populations, Tensor<double, 3> &after_collision_populations, const std::vector<BoundaryPoint<2>> &open_boundary_nodes)
    {

        Eigen::Index j, i;

        auto n_rows = populations.dimensions()[0];
        auto n_cols = populations.dimensions()[1];

        for(const auto& open_boundary : open_boundary_nodes)
        {
            i = open_boundary.coords[0];
            j = open_boundary.coords[1];

            populations(i, j, 0) = after_collision_populations(i,j,0);
            
            // If the node is not on the left wall we can propagate to the left
            if(j != 0)
            {
                populations(i, j-1, 3) = after_collision_populations(i, j, 3);
            }
            // If the node is not on the bottom wall we can propagate downwards
            if(i != n_rows-1)
            {
                populations(i+1, j, 4) = after_collision_populations(i, j, 4);
            }
            // If the node is not on the right wall we can propagate to the right
            if(j != n_cols-1)
            {
                populations(i, j+1, 1) = after_collision_populations(i, j, 1);
            }
            // If the node is not on the top wall we can propagate upwards
            if(i != 0)
            {
                populations(i-1, j, 2) = after_collision_populations(i, j, 2);
            }
            // If the node is not on the bottom left corner we can propagate to the bottom left
            if(i != n_rows-1 && j != 0)
            {
                populations(i+1, j-1, 7) = after_collision_populations(i, j, 7);
            }
            // If the node is not on the bottom right corner we can propagate to the bottom right
            if(i != n_rows-1 && j != n_cols-1)
            {
                populations(i+1, j+1, 8) = after_collision_populations(i, j, 8);
            }
            // If the node is not on the top right corner we can propagate to the top right
            if(i != 0 && j != n_cols-1)
            {
                populations(i-1, j+1, 5) = after_collision_populations(i, j, 5);
            }
            // If the node is not on the top left corner we can propagate to the top left
            if(i != 0 && j != 0)
            {
                populations(i-1, j-1, 6) = after_collision_populations(i, j, 6);
            }
        }
    }

    };


    template<std::size_t dim>
    class OMPTRTCollisionPolicy : public CollisionPolicyTag, public OMPTag
    {
    private:
        /// @brief Logger of the collision policy
        Logger l;

    public:
        OMPTRTCollisionPolicy() 
        : l ("GENERIC TRT", std::cout)
        {
            l.error("GENERIC TRT COLLISION NOT IMPLEMENTED");
            assert(false);
        }
    };


    template<>
    class OMPTRTCollisionPolicy<2> : public CollisionPolicyTag, public OMPTag
    {

    public:
        static double tau_even;
        static double tau_odd;
        static double sound_velocity;

        /**
         * @brief Initializes the TRT collision policy with the
         * relaxation constants.
         * 
         * @param tau_even the relaxation constant associated to even moments
         * @param tau_odd the relaxation constant associated to odd moments
         * @param delta_t the time step of the simulation
         * 
         */
        static void initialize(const double& tau_even_, const double& tau_odd_, const double& sound_velocity_)
        {
            tau_even = tau_even_;
            tau_odd = tau_odd_;
            sound_velocity = sound_velocity_;

            std::cout << "Inizialization done:" << std::endl;
            std::cout << " Tau+ : " << tau_even << std::endl;
            std::cout << " Tau- : " << tau_odd << std::endl;
            std::cout << " Sound velocity: " << sound_velocity << std::endl;
        }

        /**
         * @brief Computes the kinematic shear viscosity of the fluid
         * using the parameters of the collision model.
         * 
         * @return the kinematic viscosity of the fluid
         */
        static double compute_kinematic_viscosity()
        {
            return (tau_even - 0.5) * sound_velocity * sound_velocity;
        }

        /**
         * @brief Computes the magic parameter of the model, which shows
         * how much the relaxation constants control stability and accuracy
         * 
         * @note Lambda = 1/12 cancels the third order spatial error
         * @note Lambda =  1/6  cancels the fourth order spatial error
         * @note Lambda =  1/4  provides stable simulations
         * 
         * @return the magic parameter of the model, Lambda 
         */
        static double compute_magic_parameter()
        {
            std::cout <<  "Magic Parameter: " << (tau_even - 0.5)*(tau_odd- 0.5) << std::endl;
            return (tau_even - 0.5)*(tau_odd- 0.5);
        }

        /**
         * @brief Sets the value of te free paramer for odd velocities by enforcing a magic parameer value,
         * which is used to control the stability and accuracy of the simulation, and is set to 1/4 by default.
         * 
         * @param magic_parameter magic parameter target
         */
        static void enforce_magic_parameter(const double& magic_parameter = 1.0/4.0)
        {
            tau_odd = (magic_parameter)*(1.0/(tau_even-0.5))+0.5;
            std::cout << "Magic Parameter forced to " << magic_parameter << std::endl;
            std::cout << " tau+ : " << tau_even << std::endl;
            std::cout << " tau- : " << tau_odd << std::endl;
        }

        /**
         * @brief Computes the collision of populations of fluid nodes
         * 
         * @param populations of the node
         * @param equilibrium_populations of the node
         * @param after_collision_populations of the node to be streamed
         * @param fluid_nodes vector containing coordinates of all the fluid nodes of the lattice
         * @param global_rho tensor containing all densities
         * @param global_u tensor containing all velocities 
         */
        static void collide(const Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u, const double time_step)
        {

            Eigen::Index j, i;

            #pragma omp parallel for private(j, i)
            for(const auto& fluid_node : fluid_nodes)
            {
                i = fluid_node.coords[0];
                j = fluid_node.coords[1];
                after_collision_populations(i,j,0) = populations(i,j,0) - 
                    (1.0/tau_even) * (populations(i,j,0) - equilibrium_populations(i,j,0));
                
                after_collision_populations(i,j,1) = populations(i,j,1) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,1) + populations(i,j,3) - equilibrium_populations(i,j,1) - equilibrium_populations(i,j,3)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,1) - populations(i,j,3) - equilibrium_populations(i,j,1) + equilibrium_populations(i,j,3));
                
                after_collision_populations(i,j,2) = populations(i,j,2) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,2) + populations(i,j,4) - equilibrium_populations(i,j,2) - equilibrium_populations(i,j,4)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,2) - populations(i,j,4) - equilibrium_populations(i,j,2) + equilibrium_populations(i,j,4));
                
                after_collision_populations(i,j,3) = populations(i,j,3) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,3) + populations(i,j,1) - equilibrium_populations(i,j,3) - equilibrium_populations(i,j,1)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,3) - populations(i,j,1) - equilibrium_populations(i,j,3) + equilibrium_populations(i,j,1));
                
                after_collision_populations(i,j,4) = populations(i,j,4) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,4) + populations(i,j,2) - equilibrium_populations(i,j,4) - equilibrium_populations(i,j,2)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,4) - populations(i,j,2) - equilibrium_populations(i,j,4) + equilibrium_populations(i,j,2));
                
                after_collision_populations(i,j,5) = populations(i,j,5) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,5) + populations(i,j,7) - equilibrium_populations(i,j,5) - equilibrium_populations(i,j,7)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,5) - populations(i,j,7) - equilibrium_populations(i,j,5) + equilibrium_populations(i,j,7));
                
                after_collision_populations(i,j,6) = populations(i,j,6) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,6) + populations(i,j,8) - equilibrium_populations(i,j,6) - equilibrium_populations(i,j,8)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,6) - populations(i,j,8) - equilibrium_populations(i,j,6) + equilibrium_populations(i,j,8));
                
                after_collision_populations(i,j,7) = populations(i,j,7) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,7) + populations(i,j,5) - equilibrium_populations(i,j,7) - equilibrium_populations(i,j,5)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,7) - populations(i,j,5) - equilibrium_populations(i,j,7) + equilibrium_populations(i,j,5));
                
                after_collision_populations(i,j,8) = populations(i,j,8) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,8) + populations(i,j,6) - equilibrium_populations(i,j,8) - equilibrium_populations(i,j,6)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,8) - populations(i,j,6) - equilibrium_populations(i,j,8) + equilibrium_populations(i,j,6));
                
            }
        }

        /**
         * @brief Performs the streaming from after_collision_populations to populations for every fluid node
         * 
         * @param populations of fluid nodes to be updated
         * @param after_collision_populations of fluid nodes to be streamed
         * @param fluid_nodes vector containing coordinates of all the fluid nodes of the lattice
         */
        static void stream(Tensor<double, 3> &populations, Tensor<double, 3> &after_collision_populations, const std::vector<Point<2>> &fluid_nodes)
        {

            Eigen::Index j, i;

            #pragma omp parallel for private (j, i)
            for(const auto& fluid_node : fluid_nodes)
            {
                i = fluid_node.coords[0];
                j = fluid_node.coords[1];
                
                populations(i, j, 0) = after_collision_populations(i,j,0);
                populations(i, j+1, 1) = after_collision_populations(i,j,1);
                populations(i-1, j, 2) = after_collision_populations(i,j,2);
                populations(i, j-1, 3) = after_collision_populations(i,j,3);
                populations(i+1, j, 4) = after_collision_populations(i,j,4);
                populations(i-1, j+1, 5) = after_collision_populations(i,j,5);
                populations(i-1, j-1, 6) = after_collision_populations(i,j,6);
                populations(i+1, j-1, 7) = after_collision_populations(i,j,7);
                populations(i+1, j+1, 8) = after_collision_populations(i,j,8);
            }
        }

    /**
     * @brief Computes the collision of populations at the open boundaries
     * 
     * @param populations of the node
     * @param equilibrium_populations of the node
     * @param after_collision_populations of the node to be streamed
     * @param open_boundary_nodes vector containing coordinates of all the border nodes of the lattice
     * @param global_rho tensor containing all densities
     * @param global_u tensor containing all velocities 
     */
    static void collide_open_boundary(const Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, const std::vector<BoundaryPoint<2>> &open_boundary_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u, const double time_step)
    {

        Eigen::Index j, i;

        #pragma omp parallel for private(j, i)
        for(const auto& open_boundary : open_boundary_nodes)
        {
            i = open_boundary.coords[0];
            j = open_boundary.coords[1];
            after_collision_populations(i,j,0) = populations(i,j,0) - 
                (1.0/tau_even) * (populations(i,j,0) - equilibrium_populations(i,j,0));
            
            after_collision_populations(i,j,1) = populations(i,j,1) - 
                0.5 * (1.0/tau_even) * (populations(i,j,1) + populations(i,j,3) - equilibrium_populations(i,j,1) - equilibrium_populations(i,j,3)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,1) - populations(i,j,3) - equilibrium_populations(i,j,1) + equilibrium_populations(i,j,3));
            
            after_collision_populations(i,j,2) = populations(i,j,2) - 
                0.5 * (1.0/tau_even) * (populations(i,j,2) + populations(i,j,4) - equilibrium_populations(i,j,2) - equilibrium_populations(i,j,4)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,2) - populations(i,j,4) - equilibrium_populations(i,j,2) + equilibrium_populations(i,j,4));
            
            after_collision_populations(i,j,3) = populations(i,j,3) - 
                0.5 * (1.0/tau_even) * (populations(i,j,3) + populations(i,j,1) - equilibrium_populations(i,j,3) - equilibrium_populations(i,j,1)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,3) - populations(i,j,1) - equilibrium_populations(i,j,3) + equilibrium_populations(i,j,1));
            
            after_collision_populations(i,j,4) = populations(i,j,4) - 
                0.5 * (1.0/tau_even) * (populations(i,j,4) + populations(i,j,2) - equilibrium_populations(i,j,4) - equilibrium_populations(i,j,2)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,4) - populations(i,j,2) - equilibrium_populations(i,j,4) + equilibrium_populations(i,j,2));
            
            after_collision_populations(i,j,5) = populations(i,j,5) - 
                0.5 * (1.0/tau_even) * (populations(i,j,5) + populations(i,j,7) - equilibrium_populations(i,j,5) - equilibrium_populations(i,j,7)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,5) - populations(i,j,7) - equilibrium_populations(i,j,5) + equilibrium_populations(i,j,7));
            
            after_collision_populations(i,j,6) = populations(i,j,6) - 
                0.5 * (1.0/tau_even) * (populations(i,j,6) + populations(i,j,8) - equilibrium_populations(i,j,6) - equilibrium_populations(i,j,8)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,6) - populations(i,j,8) - equilibrium_populations(i,j,6) + equilibrium_populations(i,j,8));
            
            after_collision_populations(i,j,7) = populations(i,j,7) - 
                0.5 * (1.0/tau_even) * (populations(i,j,7) + populations(i,j,5) - equilibrium_populations(i,j,7) - equilibrium_populations(i,j,5)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,7) - populations(i,j,5) - equilibrium_populations(i,j,7) + equilibrium_populations(i,j,5));
            
            after_collision_populations(i,j,8) = populations(i,j,8) - 
                0.5 * (1.0/tau_even) * (populations(i,j,8) + populations(i,j,6) - equilibrium_populations(i,j,8) - equilibrium_populations(i,j,6)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,8) - populations(i,j,6) - equilibrium_populations(i,j,8) + equilibrium_populations(i,j,6));
            
        }
    }

    /**
     * @brief Performs the streaming from after_collision_populations to populations for every open boundary node
     * 
     * @param populations of border nodes to be updated
     * @param after_collision_populations of border nodes to be streamed
     * @param open_boundary_nodes vector containing coordinates of all the border nodes of the lattice
     */
    static void stream_open_boundary(Tensor<double, 3> &populations, Tensor<double, 3> &after_collision_populations, const std::vector<BoundaryPoint<2>> &open_boundary_nodes)
    {

        Eigen::Index j, i;

        auto n_rows = populations.dimensions()[0];
        auto n_cols = populations.dimensions()[1];

        #pragma omp parallel for private(j, i)
        for(const auto& open_boundary : open_boundary_nodes)
        {
            i = open_boundary.coords[0];
            j = open_boundary.coords[1];

            populations(i, j, 0) = after_collision_populations(i,j,0);
            
            // If the node is not on the left wall we can propagate to the left
            if(j != 0)
            {
                populations(i, j-1, 3) = after_collision_populations(i, j, 3);
            }
            // If the node is not on the bottom wall we can propagate downwards
            if(i != n_rows-1)
            {
                populations(i+1, j, 4) = after_collision_populations(i, j, 4);
            }
            // If the node is not on the right wall we can propagate to the right
            if(j != n_cols-1)
            {
                populations(i, j+1, 1) = after_collision_populations(i, j, 1);
            }
            // If the node is not on the top wall we can propagate upwards
            if(i != 0)
            {
                populations(i-1, j, 2) = after_collision_populations(i, j, 2);
            }
            // If the node is not on the bottom left corner we can propagate to the bottom left
            if(i != n_rows-1 && j != 0)
            {
                populations(i+1, j-1, 7) = after_collision_populations(i, j, 7);
            }
            // If the node is not on the bottom right corner we can propagate to the bottom right
            if(i != n_rows-1 && j != n_cols-1)
            {
                populations(i+1, j+1, 8) = after_collision_populations(i, j, 8);
            }
            // If the node is not on the top right corner we can propagate to the top right
            if(i != 0 && j != n_cols-1)
            {
                populations(i-1, j+1, 5) = after_collision_populations(i, j, 5);
            }
            // If the node is not on the top left corner we can propagate to the top left
            if(i != 0 && j != 0)
            {
                populations(i-1, j-1, 6) = after_collision_populations(i, j, 6);
            }
        }
    }

    };


        /**
     * @brief Most generic version of the TRT collision operator.
     * 
     * @tparam dim the total number of spatial dimensions.
     */
    template<std::size_t dim>
    class STDExecTRTCollisionPolicy : public CollisionPolicyTag, public STDExecTag
    {
    private:
        /// @brief Logger of the collision policy
        Logger l;

    public:
        STDExecTRTCollisionPolicy() 
        : l ("GENERIC TRT", std::cout)
        {
            l.error("GENERIC TRT COLLISION NOT IMPLEMENTED");
            assert(false);
        }
    };


    /**
     * @brief Specialization of the TRT collision operator for the 
     * use-case of two spatial dimensions.
     */
    template<>
    class STDExecTRTCollisionPolicy<2> : public CollisionPolicyTag, public STDExecTag
    {

    public:
        static double tau_even;
        static double tau_odd;
        static double sound_velocity;

        /**
         * @brief Initializes the TRT collision policy with the
         * relaxation constants.
         * 
         * @param tau_even the relaxation constant associated to even moments
         * @param tau_odd the relaxation constant associated to odd moments
         * @param delta_t the time step of the simulation
         * 
         */
        static void initialize(const double& tau_even_, const double& tau_odd_, const double& sound_velocity_)
        {
            tau_even = tau_even_;
            tau_odd = tau_odd_;
            sound_velocity = sound_velocity_;

            std::cout << "Inizialization done:" << std::endl;
            std::cout << " Tau+ : " << tau_even << std::endl;
            std::cout << " Tau- : " << tau_odd << std::endl;
            std::cout << " Sound velocity: " << sound_velocity << std::endl;
        }

        /**
         * @brief Computes the kinematic shear viscosity of the fluid
         * using the parameters of the collision model.
         * 
         * @return the kinematic viscosity of the fluid
         */
        static double compute_kinematic_viscosity()
        {
            return (tau_even - 0.5) * sound_velocity * sound_velocity;
        }

        /**
         * @brief Computes the magic parameter of the model, which shows
         * how much the relaxation constants control stability and accuracy
         * 
         * @note Lambda = 1/12 cancels the third order spatial error
         * @note Lambda =  1/6  cancels the fourth order spatial error
         * @note Lambda =  1/4  provides stable simulations
         * 
         * @return the magic parameter of the model, Lambda 
         */
        static double compute_magic_parameter()
        {
            std::cout <<  "Magic Parameter: " << (tau_even - 0.5)*(tau_odd- 0.5) << std::endl;
            return (tau_even - 0.5)*(tau_odd- 0.5);
        }

        /**
         * @brief Sets the value of te free paramer for odd velocities by enforcing a magic parameer value,
         * which is used to control the stability and accuracy of the simulation, and is set to 1/4 by default.
         * 
         * @param magic_parameter magic parameter target
         */
        static void enforce_magic_parameter(const double& magic_parameter = 1.0/4.0)
        {
            tau_odd = (magic_parameter)*(1.0/(tau_even-0.5))+0.5;
            std::cout << "Magic Parameter forced to " << magic_parameter << std::endl;
            std::cout << " tau+ : " << tau_even << std::endl;
            std::cout << " tau- : " << tau_odd << std::endl;
        }

        /**
         * @brief Computes the collision of populations of fluid nodes
         * 
         * @param populations of the node
         * @param equilibrium_populations of the node
         * @param after_collision_populations of the node to be streamed
         * @param fluid_nodes vector containing coordinates of all the fluid nodes of the lattice
         * @param global_rho tensor containing all densities
         * @param global_u tensor containing all velocities 
         */
        static void collide(const Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u, const double time_step)
        {

            Eigen::Index j, i;

            //for(const auto& fluid_node : fluid_nodes)
            std::for_each(std::execution::par,fluid_nodes.begin(),fluid_nodes.end(), [&](const auto fluid_node){
                i = fluid_node.coords[0];
                j = fluid_node.coords[1];
                after_collision_populations(i,j,0) = populations(i,j,0) - 
                    (1.0/tau_even) * (populations(i,j,0) - equilibrium_populations(i,j,0));
                
                after_collision_populations(i,j,1) = populations(i,j,1) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,1) + populations(i,j,3) - equilibrium_populations(i,j,1) - equilibrium_populations(i,j,3)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,1) - populations(i,j,3) - equilibrium_populations(i,j,1) + equilibrium_populations(i,j,3));
                
                after_collision_populations(i,j,2) = populations(i,j,2) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,2) + populations(i,j,4) - equilibrium_populations(i,j,2) - equilibrium_populations(i,j,4)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,2) - populations(i,j,4) - equilibrium_populations(i,j,2) + equilibrium_populations(i,j,4));
                
                after_collision_populations(i,j,3) = populations(i,j,3) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,3) + populations(i,j,1) - equilibrium_populations(i,j,3) - equilibrium_populations(i,j,1)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,3) - populations(i,j,1) - equilibrium_populations(i,j,3) + equilibrium_populations(i,j,1));
                
                after_collision_populations(i,j,4) = populations(i,j,4) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,4) + populations(i,j,2) - equilibrium_populations(i,j,4) - equilibrium_populations(i,j,2)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,4) - populations(i,j,2) - equilibrium_populations(i,j,4) + equilibrium_populations(i,j,2));
                
                after_collision_populations(i,j,5) = populations(i,j,5) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,5) + populations(i,j,7) - equilibrium_populations(i,j,5) - equilibrium_populations(i,j,7)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,5) - populations(i,j,7) - equilibrium_populations(i,j,5) + equilibrium_populations(i,j,7));
                
                after_collision_populations(i,j,6) = populations(i,j,6) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,6) + populations(i,j,8) - equilibrium_populations(i,j,6) - equilibrium_populations(i,j,8)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,6) - populations(i,j,8) - equilibrium_populations(i,j,6) + equilibrium_populations(i,j,8));
                
                after_collision_populations(i,j,7) = populations(i,j,7) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,7) + populations(i,j,5) - equilibrium_populations(i,j,7) - equilibrium_populations(i,j,5)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,7) - populations(i,j,5) - equilibrium_populations(i,j,7) + equilibrium_populations(i,j,5));
                
                after_collision_populations(i,j,8) = populations(i,j,8) - 
                    0.5 * (1.0/tau_even) * (populations(i,j,8) + populations(i,j,6) - equilibrium_populations(i,j,8) - equilibrium_populations(i,j,6)) -
                    0.5 * (1.0/tau_odd) * (populations(i,j,8) - populations(i,j,6) - equilibrium_populations(i,j,8) + equilibrium_populations(i,j,6));
                
            });
        }

        /**
         * @brief Performs the streaming from after_collision_populations to populations for every fluid node
         * 
         * @param populations of fluid nodes to be updated
         * @param after_collision_populations of fluid nodes to be streamed
         * @param fluid_nodes vector containing coordinates of all the fluid nodes of the lattice
         */
        static void stream(Tensor<double, 3> &populations, Tensor<double, 3> &after_collision_populations, const std::vector<Point<2>> &fluid_nodes)
        {

            Eigen::Index j, i;

            
            //for(const auto& fluid_node : fluid_nodes)
            std::for_each(std::execution::par,fluid_nodes.begin(),fluid_nodes.end(), [&](const auto fluid_node)
            {
                i = fluid_node.coords[0];
                j = fluid_node.coords[1];
                
                populations(i, j, 0) = after_collision_populations(i,j,0);
                populations(i, j+1, 1) = after_collision_populations(i,j,1);
                populations(i-1, j, 2) = after_collision_populations(i,j,2);
                populations(i, j-1, 3) = after_collision_populations(i,j,3);
                populations(i+1, j, 4) = after_collision_populations(i,j,4);
                populations(i-1, j+1, 5) = after_collision_populations(i,j,5);
                populations(i-1, j-1, 6) = after_collision_populations(i,j,6);
                populations(i+1, j-1, 7) = after_collision_populations(i,j,7);
                populations(i+1, j+1, 8) = after_collision_populations(i,j,8);
            });
        }

    /**
     * @brief Computes the collision of populations at the open boundaries
     * 
     * @param populations of the node
     * @param equilibrium_populations of the node
     * @param after_collision_populations of the node to be streamed
     * @param open_boundary_nodes vector containing coordinates of all the border nodes of the lattice
     * @param global_rho tensor containing all densities
     * @param global_u tensor containing all velocities 
     */
    static void collide_open_boundary(const Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &after_collision_populations, const std::vector<BoundaryPoint<2>> &open_boundary_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u, const double time_step)
    {

        Eigen::Index j, i;

        //for(const auto& open_boundary : open_boundary_nodes)
        std::for_each(std::execution::par,open_boundary_nodes.begin(),open_boundary_nodes.end(), [&](const auto open_boundary)
        {
            i = open_boundary.coords[0];
            j = open_boundary.coords[1];
            after_collision_populations(i,j,0) = populations(i,j,0) - 
                (1.0/tau_even) * (populations(i,j,0) - equilibrium_populations(i,j,0));
            
            after_collision_populations(i,j,1) = populations(i,j,1) - 
                0.5 * (1.0/tau_even) * (populations(i,j,1) + populations(i,j,3) - equilibrium_populations(i,j,1) - equilibrium_populations(i,j,3)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,1) - populations(i,j,3) - equilibrium_populations(i,j,1) + equilibrium_populations(i,j,3));
            
            after_collision_populations(i,j,2) = populations(i,j,2) - 
                0.5 * (1.0/tau_even) * (populations(i,j,2) + populations(i,j,4) - equilibrium_populations(i,j,2) - equilibrium_populations(i,j,4)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,2) - populations(i,j,4) - equilibrium_populations(i,j,2) + equilibrium_populations(i,j,4));
            
            after_collision_populations(i,j,3) = populations(i,j,3) - 
                0.5 * (1.0/tau_even) * (populations(i,j,3) + populations(i,j,1) - equilibrium_populations(i,j,3) - equilibrium_populations(i,j,1)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,3) - populations(i,j,1) - equilibrium_populations(i,j,3) + equilibrium_populations(i,j,1));
            
            after_collision_populations(i,j,4) = populations(i,j,4) - 
                0.5 * (1.0/tau_even) * (populations(i,j,4) + populations(i,j,2) - equilibrium_populations(i,j,4) - equilibrium_populations(i,j,2)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,4) - populations(i,j,2) - equilibrium_populations(i,j,4) + equilibrium_populations(i,j,2));
            
            after_collision_populations(i,j,5) = populations(i,j,5) - 
                0.5 * (1.0/tau_even) * (populations(i,j,5) + populations(i,j,7) - equilibrium_populations(i,j,5) - equilibrium_populations(i,j,7)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,5) - populations(i,j,7) - equilibrium_populations(i,j,5) + equilibrium_populations(i,j,7));
            
            after_collision_populations(i,j,6) = populations(i,j,6) - 
                0.5 * (1.0/tau_even) * (populations(i,j,6) + populations(i,j,8) - equilibrium_populations(i,j,6) - equilibrium_populations(i,j,8)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,6) - populations(i,j,8) - equilibrium_populations(i,j,6) + equilibrium_populations(i,j,8));
            
            after_collision_populations(i,j,7) = populations(i,j,7) - 
                0.5 * (1.0/tau_even) * (populations(i,j,7) + populations(i,j,5) - equilibrium_populations(i,j,7) - equilibrium_populations(i,j,5)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,7) - populations(i,j,5) - equilibrium_populations(i,j,7) + equilibrium_populations(i,j,5));
            
            after_collision_populations(i,j,8) = populations(i,j,8) - 
                0.5 * (1.0/tau_even) * (populations(i,j,8) + populations(i,j,6) - equilibrium_populations(i,j,8) - equilibrium_populations(i,j,6)) -
                0.5 * (1.0/tau_odd) * (populations(i,j,8) - populations(i,j,6) - equilibrium_populations(i,j,8) + equilibrium_populations(i,j,6));
            
        });
    }

    /**
     * @brief Performs the streaming from after_collision_populations to populations for every open boundary node
     * 
     * @param populations of border nodes to be updated
     * @param after_collision_populations of border nodes to be streamed
     * @param open_boundary_nodes vector containing coordinates of all the border nodes of the lattice
     */
    static void stream_open_boundary(Tensor<double, 3> &populations, Tensor<double, 3> &after_collision_populations, const std::vector<BoundaryPoint<2>> &open_boundary_nodes)
    {

        Eigen::Index j, i;

        auto n_rows = populations.dimensions()[0];
        auto n_cols = populations.dimensions()[1];

        std::for_each(std::execution::par,open_boundary_nodes.begin(),open_boundary_nodes.end(), [&](const auto open_boundary)
        {
            i = open_boundary.coords[0];
            j = open_boundary.coords[1];

            populations(i, j, 0) = after_collision_populations(i,j,0);
            
            // If the node is not on the left wall we can propagate to the left
            if(j != 0)
            {
                populations(i, j-1, 3) = after_collision_populations(i, j, 3);
            }
            // If the node is not on the bottom wall we can propagate downwards
            if(i != n_rows-1)
            {
                populations(i+1, j, 4) = after_collision_populations(i, j, 4);
            }
            // If the node is not on the right wall we can propagate to the right
            if(j != n_cols-1)
            {
                populations(i, j+1, 1) = after_collision_populations(i, j, 1);
            }
            // If the node is not on the top wall we can propagate upwards
            if(i != 0)
            {
                populations(i-1, j, 2) = after_collision_populations(i, j, 2);
            }
            // If the node is not on the bottom left corner we can propagate to the bottom left
            if(i != n_rows-1 && j != 0)
            {
                populations(i+1, j-1, 7) = after_collision_populations(i, j, 7);
            }
            // If the node is not on the bottom right corner we can propagate to the bottom right
            if(i != n_rows-1 && j != n_cols-1)
            {
                populations(i+1, j+1, 8) = after_collision_populations(i, j, 8);
            }
            // If the node is not on the top right corner we can propagate to the top right
            if(i != 0 && j != n_cols-1)
            {
                populations(i-1, j+1, 5) = after_collision_populations(i, j, 5);
            }
            // If the node is not on the top left corner we can propagate to the top left
            if(i != 0 && j != 0)
            {
                populations(i-1, j-1, 6) = after_collision_populations(i, j, 6);
            }
        });
    }

    };

    // initialization of the relaxation constants in the 2-D TRT collision operator.
    double TRTCollisionPolicy<2>::tau_even = 0.0;
    double TRTCollisionPolicy<2>::tau_odd = 0.0;
    double TRTCollisionPolicy<2>::sound_velocity = 0.0;


    // initialization of the relaxation constants in the 2-D TRT collision operator(for OMP class).
    double OMPTRTCollisionPolicy<2>::tau_even = 0.0;
    double OMPTRTCollisionPolicy<2>::tau_odd = 0.0;
    double OMPTRTCollisionPolicy<2>::sound_velocity = 0.0;

    // initialization of the relaxation constants in the 2-D TRT collision operator(for STDExec class).
    double STDExecTRTCollisionPolicy<2>::tau_even = 0.0;
    double STDExecTRTCollisionPolicy<2>::tau_odd = 0.0;
    double STDExecTRTCollisionPolicy<2>::sound_velocity = 0.0;
}; // namespace llalbm::core::collisions

#endif // LLALBM_TRT_COLLISION_HPP