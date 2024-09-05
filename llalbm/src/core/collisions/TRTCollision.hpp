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
#include <openacc.h>
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


        /**
     * @brief Most generic version of the TRT collision operator.
     * 
     * @tparam dim the total number of spatial dimensions.
     */
    template<std::size_t dim>
    class OpenACCTRTCollisionPolicy : public CollisionPolicyTag, public OpenACCTag
    {
    private:
        /// @brief Logger of the collision policy
        Logger l;

    public:
        OpenACCTRTCollisionPolicy() 
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
    class OpenACCTRTCollisionPolicy<2> : public CollisionPolicyTag, public OpenACCTag
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
            Eigen::Index i;
            Eigen::Index j;

            auto n_rows = populations.dimensions()[0];
            auto n_cols = populations.dimensions()[1];
            auto n_dirs = populations.dimensions()[2]; 

            double* populations_buffer = new double[n_rows * n_cols * n_dirs];
            double* equilibrium_populations_buffer = new double[n_rows * n_cols * n_dirs];
            double* after_collision_populations_buffer = new double[n_rows * n_cols * n_dirs];

            for (size_t i = 0; i < n_rows; ++i) {
                for (size_t j = 0; j < n_cols; ++j) {
                    for (size_t d = 0; d < n_dirs; ++d) {
                        populations_buffer[i * n_cols * n_dirs + j * n_dirs + d] = populations(i, j, d);
                        equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + d] = equilibrium_populations(i, j, d);
                        after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + d] = after_collision_populations(i, j, d);
                    }
                }
            }

            #pragma acc data copy(after_collision_populations_buffer[n_rows * n_cols * n_dirs]) copyin(equilibrium_populations_buffer[n_rows * n_cols * n_dirs], populations_buffer[n_rows * n_cols * n_dirs])
            #pragma acc cache(populations_buffer[n_rows * n_cols * n_dirs], equilibrium_populations_buffer[n_rows * n_cols * n_dirs], after_collision_populations_buffer[n_rows * n_cols * n_dirs])
            {
                #pragma acc parallel loop
                for (size_t fnode = 0; fnode < fluid_nodes.size(); ++fnode) {
                    i = fluid_nodes[fnode].coords[0];
                    j = fluid_nodes[fnode].coords[1];

                    after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 0] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 0] - 
                        (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 0] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 0]);

                    after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] - 
                        0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3]) -
                        0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3]);

                    after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] - 
                        0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4]) -
                        0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4]);

                    after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] - 
                        0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1]) -
                        0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1]);

                    after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] - 
                        0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2]) -
                        0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2]);

                    after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] - 
                        0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7]) -
                        0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7]);

                    after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] - 
                        0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8]) -
                        0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8]);

                    after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] - 
                        0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5]) -
                        0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5]);

                    after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] - 
                        0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6]) -
                        0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6]);
                }

                for (size_t i = 0; i < n_rows; ++i) {
                    for (size_t j = 0; j < n_cols; ++j) {
                        for (size_t d = 0; d < n_dirs; ++d) {
                            after_collision_populations(i, j, d) = after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + d];
                        }
                    }
                }
            }

            delete[] populations_buffer;
            delete[] after_collision_populations_buffer;
            delete[] equilibrium_populations_buffer;


            /*Eigen::Index j, i;
            
            #pragma acc parallel loop
            //for(const auto& fluid_node : fluid_nodes)
            for(size_t fnode=0; fnode < fluid_nodes.size(); fnode++)
            {
                i = fluid_nodes[fnode].coords[0];
                j = fluid_nodes[fnode].coords[1];
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
                
            }*/
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
            auto n_rows = populations.dimensions()[0];
            auto n_cols = populations.dimensions()[1];
            auto num_directions = populations.dimensions()[2];
            double* populations_buffer = new double[n_rows * n_cols * num_directions];
            double* after_collision_populations_buffer = new double[n_rows * n_cols * num_directions];

            // Copia i dati da after_collision_populations al buffer after_collision_populations_buffer
            for (Eigen::Index i = 0; i < n_rows; ++i) {
                for (Eigen::Index j = 0; j < n_cols; ++j) {
                    for (std::size_t d = 0; d < num_directions; ++d) {
                        after_collision_populations_buffer[i * n_cols * num_directions + j * num_directions + d] = after_collision_populations(i, j, d);
                    }
                }
            }

            // Loop parallelo per aggiornare populations
            #pragma acc data copy(populations_buffer[0:n_rows*n_cols*num_directions])
            {
                #pragma acc parallel loop
                for (size_t fnode = 0; fnode < fluid_nodes.size(); fnode++) {
                    size_t i = fluid_nodes[fnode].coords[0];
                    size_t j = fluid_nodes[fnode].coords[1];

                    // Aggiorna il buffer di populations
                    populations_buffer[i * n_cols * num_directions + j * num_directions + 0] = after_collision_populations_buffer[i * n_cols * num_directions + j * num_directions + 0];
                    populations_buffer[i * n_cols * num_directions + (j + 1) * num_directions + 1] = after_collision_populations_buffer[i * n_cols * num_directions + j * num_directions + 1];
                    populations_buffer[(i - 1) * n_cols * num_directions + j * num_directions + 2] = after_collision_populations_buffer[i * n_cols * num_directions + j * num_directions + 2];
                    populations_buffer[i * n_cols * num_directions + (j - 1) * num_directions + 3] = after_collision_populations_buffer[i * n_cols * num_directions + j * num_directions + 3];
                    populations_buffer[(i + 1) * n_cols * num_directions + j * num_directions + 4] = after_collision_populations_buffer[i * n_cols * num_directions + j * num_directions + 4];
                    populations_buffer[(i - 1) * n_cols * num_directions + (j + 1) * num_directions + 5] = after_collision_populations_buffer[i * n_cols * num_directions + j * num_directions + 5];
                    populations_buffer[(i - 1) * n_cols * num_directions + (j - 1) * num_directions + 6] = after_collision_populations_buffer[i * n_cols * num_directions + j * num_directions + 6];
                    populations_buffer[(i + 1) * n_cols * num_directions + (j - 1) * num_directions + 7] = after_collision_populations_buffer[i * n_cols * num_directions + j * num_directions + 7];
                    populations_buffer[(i + 1) * n_cols * num_directions + (j + 1) * num_directions + 8] = after_collision_populations_buffer[i * n_cols * num_directions + j * num_directions + 8];
                }
            }
            // Libera la memoria allocata dinamicamente
            delete[] populations_buffer;
            delete[] after_collision_populations_buffer;

            /*#pragma acc parallel loop
            //for(const auto& fluid_node : fluid_nodes)
            for(size_t fnode=0; fnode < fluid_nodes.size(); fnode++)
            {
                i = fluid_nodes[fnode].coords[0];
                j = fluid_nodes[fnode].coords[1];
                
                populations(i, j, 0) = after_collision_populations(i,j,0);
                populations(i, j+1, 1) = after_collision_populations(i,j,1);
                populations(i-1, j, 2) = after_collision_populations(i,j,2);
                populations(i, j-1, 3) = after_collision_populations(i,j,3);
                populations(i+1, j, 4) = after_collision_populations(i,j,4);
                populations(i-1, j+1, 5) = after_collision_populations(i,j,5);
                populations(i-1, j-1, 6) = after_collision_populations(i,j,6);
                populations(i+1, j-1, 7) = after_collision_populations(i,j,7);
                populations(i+1, j+1, 8) = after_collision_populations(i,j,8);
            }*/
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

        auto n_rows = populations.dimensions()[0];
        auto n_cols = populations.dimensions()[1];
        auto n_dirs = populations.dimensions()[2]; 

        double* populations_buffer = new double[n_rows * n_cols * n_dirs];
        double* equilibrium_populations_buffer = new double[n_rows * n_cols * n_dirs];
        double* after_collision_populations_buffer = new double[n_rows * n_cols * n_dirs];

        for (size_t i = 0; i < n_rows; ++i) {
            for (size_t j = 0; j < n_cols; ++j) {
                for (size_t d = 0; d < n_dirs; ++d) {
                    populations_buffer[i * n_cols * n_dirs + j * n_dirs + d] = populations(i, j, d);
                    equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + d] = equilibrium_populations(i, j, d);
                    after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + d] = after_collision_populations(i, j, d);
                }
            }
        }

        #pragma acc data copy(after_collision_populations_buffer[n_rows * n_cols * n_dirs]) copyin(equilibrium_populations_buffer[n_rows * n_cols * n_dirs], populations_buffer[n_rows * n_cols * n_dirs])
        #pragma acc cache(populations_buffer[n_rows * n_cols * n_dirs], equilibrium_populations_buffer[n_rows * n_cols * n_dirs], after_collision_populations_buffer[n_rows * n_cols * n_dirs])
        {
            #pragma acc parallel loop
            for (size_t onode = 0; onode < open_boundary_nodes.size(); ++onode) {
                Eigen::Index i = open_boundary_nodes[onode].coords[0];
                Eigen::Index j = open_boundary_nodes[onode].coords[1];

                after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 0] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 0] - 
                    (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 0] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 0]);

                after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] - 
                    0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3]) -
                    0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3]);

                after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] - 
                    0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4]) -
                    0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4]);

                after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] - 
                    0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1]) -
                    0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1]);

                after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] - 
                    0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2]) -
                    0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2]);

                after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] - 
                    0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7]) -
                    0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7]);

                after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] - 
                    0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8]) -
                    0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8]);

                after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] - 
                    0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5]) -
                    0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5]);

                after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] = populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] - 
                    0.5 * (1.0 / tau_even) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] + populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6]) -
                    0.5 * (1.0 / tau_odd) * (populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] - populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6] - equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8] + equilibrium_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6]);
            }

            for (size_t i = 0; i < n_rows; ++i) {
                    for (size_t j = 0; j < n_cols; ++j) {
                        for (size_t d = 0; d < n_dirs; ++d) {
                            after_collision_populations(i, j, d) = after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + d];
                        }
                    }
                }
            }

        delete[] populations_buffer;
        delete[] after_collision_populations_buffer;
        delete[] equilibrium_populations_buffer;
        /*#pragma acc parallel loop
        //for(const auto& open_boundary : open_boundary_nodes)
        for(size_t onode=0; onode < open_boundary_nodes.size(); onode++)
        {
            i = open_boundary_nodes[onode].coords[0];
            j = open_boundary_nodes[onode].coords[1];
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
            
        }*/
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
        auto n_dirs = populations.dimensions()[2]; // Presumendo che le direzioni siano la terza dimensione

        
        double* populations_buffer = new double[n_rows * n_cols * n_dirs];
        double* after_collision_populations_buffer= new double[n_rows * n_cols * n_dirs];

        
        for (size_t i = 0; i < n_rows; ++i) {
            for (size_t j = 0; j < n_cols; ++j) {
                for (size_t d = 0; d < n_dirs; ++d) {
                    populations_buffer[i * n_cols * n_dirs + j * n_dirs + d] = populations(i, j, d);
                    after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + d] = after_collision_populations(i, j, d);
                }
            }
        }

        // Uso dei buffer
        #pragma acc data copy(populations_buffer[n_rows * n_cols * n_dirs], after_collision_populations_buffer[n_rows * n_cols * n_dirs])
        {
            #pragma acc parallel loop
            for (size_t onode = 0; onode < open_boundary_nodes.size(); ++onode) {
                size_t i = open_boundary_nodes[onode].coords[0];
                size_t j = open_boundary_nodes[onode].coords[1];

                populations_buffer[i * n_cols * n_dirs + j * n_dirs + 0] = after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 0];

                if (j != 0) {
                    populations_buffer[i * n_cols * n_dirs + (j - 1) * n_dirs + 3] = after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 3];
                }
                if (i != n_rows - 1) {
                    populations_buffer[(i + 1) * n_cols * n_dirs + j * n_dirs + 4] = after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 4];
                }
                if (j != n_cols - 1) {
                    populations_buffer[i * n_cols * n_dirs + (j + 1) * n_dirs + 1] = after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 1];
                }
                if (i != 0) {
                    populations_buffer[(i - 1) * n_cols * n_dirs + j * n_dirs + 2] = after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 2];
                }
                if (i != n_rows - 1 && j != 0) {
                    populations_buffer[(i + 1) * n_cols * n_dirs + (j - 1) * n_dirs + 7] = after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 7];
                }
                if (i != n_rows - 1 && j != n_cols - 1) {
                    populations_buffer[(i + 1) * n_cols * n_dirs + (j + 1) * n_dirs + 8] = after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 8];
                }
                if (i != 0 && j != n_cols - 1) {
                    populations_buffer[(i - 1) * n_cols * n_dirs + (j + 1) * n_dirs + 5] = after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 5];
                }
                if (i != 0 && j != 0) {
                    populations_buffer[(i - 1) * n_cols * n_dirs + (j - 1) * n_dirs + 6] = after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + 6];
                }
            }
        }

        for (size_t i = 0; i < n_rows; ++i) {
            for (size_t j = 0; j < n_cols; ++j) {
                for (size_t d = 0; d < n_dirs; ++d) {
                    populations(i, j, d) = populations_buffer[i * n_cols * n_dirs + j * n_dirs + d];
                    after_collision_populations(i, j, d) = after_collision_populations_buffer[i * n_cols * n_dirs + j * n_dirs + d];
                }
            }
        }
        delete[] populations_buffer;
        delete[] after_collision_populations_buffer;

        /*auto n_rows = populations.dimensions()[0];
        auto n_cols = populations.dimensions()[1];

        #pragma acc parallel loop
        //for(const auto& open_boundary : open_boundary_nodes)
        for(size_t onode=0; onode < open_boundary_nodes.size(); onode++)
        {
            i = open_boundary_nodes[onode].coords[0];
            j = open_boundary_nodes[onode].coords[1];

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
        }*/
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

    // initialization of the relaxation constants in the 2-D TRT collision operator.
    double OpenACCTRTCollisionPolicy<2>::tau_even = 0.0;
    double OpenACCTRTCollisionPolicy<2>::tau_odd = 0.0;
    double OpenACCTRTCollisionPolicy<2>::sound_velocity = 0.0;

}; // namespace llalbm::core::collisions

#endif // LLALBM_TRT_COLLISION_HPP