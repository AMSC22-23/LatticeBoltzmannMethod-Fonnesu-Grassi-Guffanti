/**
 * @file BGKCollision.hpp File containing the templates for the BGK collision policy
 * @brief The BGK collision in the Lattice-Boltzmann method is among the 
 * 
 */

#ifndef LLALBM_BGK_COLLISION_HPP
#define LLALBM_BGK_COLLISION_HPP

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
     * @brief Most generic version of the BGK collision operator.
     * 
     * @tparam dim the total number of spatial dimensions.
     */
    template<std::size_t dim>
    class BGKCollisionPolicy : public CollisionPolicyTag, public SequentialTag
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
    class BGKCollisionPolicy<2> : public CollisionPolicyTag, public SequentialTag
    {

    public:
        static double tau;
        
        /**
         * @brief Initializes the BGK collision policy with the
         * relaxation constant.
         * 
         * @param tau_ the relaxation constant
         */
        static void initialize(const double& tau_, const double& max_speed = 0.2)
        {
            tau = tau_;
            std::cout << "Inizialization done:" << std::endl;
            std::cout << " Tau: " << tau << std::endl;
            std::cout << " Max Speed: " << max_speed << std::endl;
            std::cout << "====================================================" << std::endl;
            if (tau <= 0.55)
            {
                std::cout << " The simulation might have stability issues." << std::endl;
                std::cout << " Theoretical LIMIT of tau:" << "tau > " << 0.5 + (1/8.0) * max_speed << std::endl;
            }
            std::cout << " Allowed Reynolds Number: " << max_speed * 3.0 / (tau - 0.5) << std::endl;
        }

        /**
         * @brief Inializes the value of the relaxation parameter based on the desired maximum velocity and Reynolds number
         * of the simulation. 
         * @note Quantities are expressed in lattice units.
         * 
         * @param max_velocity maximum velocity of the simulation 
         * @param reynolds Reynolds number of the simulation
         */
        static void allowed_tau(const double& max_speed, const double& reynolds, const double& speed_of_sound_squared = 1.0/3.0)
        {
            tau = 0.5 + max_speed / (speed_of_sound_squared * reynolds);
            std::cout << "Inizialization done:" << std::endl;
            std::cout << " Max Speed: " << max_speed << std::endl;
            std::cout << " Reynolds: " << reynolds << std::endl;
            std::cout << " Speed of Sound Squared: " << speed_of_sound_squared << std::endl;
            std::cout << "====================================================" << std::endl;
            std::cout << " Smallest theoretical Tau: " << tau << std::endl;
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
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double t;
            double o_m_t;

            for(const auto& fluid_node : fluid_nodes)
            {
                i = fluid_node.coords[0];
                j = fluid_node.coords[1];

                p0 = populations(i,j,0);
                p1 = populations(i,j,1);
                p2 = populations(i,j,2);
                p3 = populations(i,j,3);
                p4 = populations(i,j,4);
                p5 = populations(i,j,5);
                p6 = populations(i,j,6);
                p7 = populations(i,j,7);
                p8 = populations(i,j,8);

                t = 1.0 / tau;  // TIME STEP IS IN LATTICE UNITS PORCO ZIO
                o_m_t = 1.0 - t;

                after_collision_populations(i,j,0) = o_m_t * p0 + t * equilibrium_populations(i,j,0);
                after_collision_populations(i,j,1) = o_m_t * p1 + t * equilibrium_populations(i,j,1);
                after_collision_populations(i,j,2) = o_m_t * p2 + t * equilibrium_populations(i,j,2);
                after_collision_populations(i,j,3) = o_m_t * p3 + t * equilibrium_populations(i,j,3);
                after_collision_populations(i,j,4) = o_m_t * p4 + t * equilibrium_populations(i,j,4);
                after_collision_populations(i,j,5) = o_m_t * p5 + t * equilibrium_populations(i,j,5);
                after_collision_populations(i,j,6) = o_m_t * p6 + t * equilibrium_populations(i,j,6);
                after_collision_populations(i,j,7) = o_m_t * p7 + t * equilibrium_populations(i,j,7);
                after_collision_populations(i,j,8) = o_m_t * p8 + t * equilibrium_populations(i,j,8);
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
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double t;
            double o_m_t;

            auto n_rows = populations.dimensions()[0];
            auto n_cols = populations.dimensions()[1];
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
        double p0, p1, p2, p3, p4, p5, p6, p7, p8;
        double t;
        double o_m_t;

        for(const auto& open_boundary : open_boundary_nodes)
        {
            i = open_boundary.coords[0];
            j = open_boundary.coords[1];

            p0 = populations(i,j,0);
            p1 = populations(i,j,1);
            p2 = populations(i,j,2);
            p3 = populations(i,j,3);
            p4 = populations(i,j,4);
            p5 = populations(i,j,5);
            p6 = populations(i,j,6);
            p7 = populations(i,j,7);
            p8 = populations(i,j,8);

            t = 1.0 / tau;  // TIME STEP IS IN LATTICE UNITS PORCO ZIO
            o_m_t = 1.0 - t;

            after_collision_populations(i,j,0) = o_m_t * p0 + t * equilibrium_populations(i,j,0);
            after_collision_populations(i,j,1) = o_m_t * p1 + t * equilibrium_populations(i,j,1);
            after_collision_populations(i,j,2) = o_m_t * p2 + t * equilibrium_populations(i,j,2);
            after_collision_populations(i,j,3) = o_m_t * p3 + t * equilibrium_populations(i,j,3);
            after_collision_populations(i,j,4) = o_m_t * p4 + t * equilibrium_populations(i,j,4);
            after_collision_populations(i,j,5) = o_m_t * p5 + t * equilibrium_populations(i,j,5);
            after_collision_populations(i,j,6) = o_m_t * p6 + t * equilibrium_populations(i,j,6);
            after_collision_populations(i,j,7) = o_m_t * p7 + t * equilibrium_populations(i,j,7);
            after_collision_populations(i,j,8) = o_m_t * p8 + t * equilibrium_populations(i,j,8);
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
        double p0, p1, p2, p3, p4, p5, p6, p7, p8;
        double t;
        double o_m_t;

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


/**
     * @brief OMP version of the BGK collision operator.
     * 
     * @tparam dim the total number of spatial dimensions.
     */
    template<std::size_t dim>
    class OMPBGKCollisionPolicy : public CollisionPolicyTag, public OMPTag
    {
    private:
        /// @brief Logger of the collision policy
        Logger l;

    public:
        OMPBGKCollisionPolicy() 
        : l ("OMP BGK", std::cout)
        {
            l.error("OMP BGK COLLISION NOT IMPLEMENTED");
            assert(false);
        }
    };


    /**
     * @brief Specialization of the omp BGK collision operator for the 
     * use-case of two spatial dimensions.
     */
    template<>
    class OMPBGKCollisionPolicy<2> : public CollisionPolicyTag, public OMPTag
    {

    public:
        static double tau;
        
        /**
         * @brief Initializes the BGK collision policy with the
         * relaxation constant.
         * 
         * @param tau_ the relaxation constant
         */
        static void initialize(const double& tau_, const double& max_speed = 0.2)
        {
            tau = tau_;
            std::cout << "Inizialization done:" << std::endl;
            std::cout << " Tau: " << tau << std::endl;
            std::cout << " Max Speed: " << max_speed << std::endl;
            std::cout << "====================================================" << std::endl;
            if (tau <= 0.55)
            {
                std::cout << " The simulation might have stability issues." << std::endl;
                std::cout << " Theoretical LIMIT of tau:" << "tau > " << 0.5 + (1/8.0) * max_speed << std::endl;
            }
            std::cout << " Allowed Reynolds Number: " << max_speed * 3.0 / (tau - 0.5) << std::endl;
        }

        /**
         * @brief Inializes the value of the relaxation parameter based on the desired maximum velocity and Reynolds number
         * of the simulation. 
         * @note Quantities are expressed in lattice units.
         * 
         * @param max_velocity maximum velocity of the simulation 
         * @param reynolds Reynolds number of the simulation
         */
        static void allowed_tau(const double& max_speed, const double& reynolds, const double& speed_of_sound_squared = 1.0/3.0)
        {
            tau = 0.5 + max_speed / (speed_of_sound_squared * reynolds);
            std::cout << "Inizialization done:" << std::endl;
            std::cout << " Max Speed: " << max_speed << std::endl;
            std::cout << " Reynolds: " << reynolds << std::endl;
            std::cout << " Speed of Sound Squared: " << speed_of_sound_squared << std::endl;
            std::cout << "====================================================" << std::endl;
            std::cout << " Smallest theoretical Tau: " << tau << std::endl;
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
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double t;
            double o_m_t;

            #pragma omp parallel for
            for(const auto& fluid_node : fluid_nodes)
            {
                i = fluid_node.coords[0];
                j = fluid_node.coords[1];

                p0 = populations(i,j,0);
                p1 = populations(i,j,1);
                p2 = populations(i,j,2);
                p3 = populations(i,j,3);
                p4 = populations(i,j,4);
                p5 = populations(i,j,5);
                p6 = populations(i,j,6);
                p7 = populations(i,j,7);
                p8 = populations(i,j,8);

                t = 1.0 / tau;  // TIME STEP IS IN LATTICE UNITS PORCO ZIO
                o_m_t = 1.0 - t;

                after_collision_populations(i,j,0) = o_m_t * p0 + t * equilibrium_populations(i,j,0);
                after_collision_populations(i,j,1) = o_m_t * p1 + t * equilibrium_populations(i,j,1);
                after_collision_populations(i,j,2) = o_m_t * p2 + t * equilibrium_populations(i,j,2);
                after_collision_populations(i,j,3) = o_m_t * p3 + t * equilibrium_populations(i,j,3);
                after_collision_populations(i,j,4) = o_m_t * p4 + t * equilibrium_populations(i,j,4);
                after_collision_populations(i,j,5) = o_m_t * p5 + t * equilibrium_populations(i,j,5);
                after_collision_populations(i,j,6) = o_m_t * p6 + t * equilibrium_populations(i,j,6);
                after_collision_populations(i,j,7) = o_m_t * p7 + t * equilibrium_populations(i,j,7);
                after_collision_populations(i,j,8) = o_m_t * p8 + t * equilibrium_populations(i,j,8);
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
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double t;
            double o_m_t;

            auto n_rows = populations.dimensions()[0];
            auto n_cols = populations.dimensions()[1];

            #pragma omp parallel for
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
        double p0, p1, p2, p3, p4, p5, p6, p7, p8;
        double t;
        double o_m_t;

        #pragma omp parallel for
        for(const auto& open_boundary : open_boundary_nodes)
        {
            i = open_boundary.coords[0];
            j = open_boundary.coords[1];

            p0 = populations(i,j,0);
            p1 = populations(i,j,1);
            p2 = populations(i,j,2);
            p3 = populations(i,j,3);
            p4 = populations(i,j,4);
            p5 = populations(i,j,5);
            p6 = populations(i,j,6);
            p7 = populations(i,j,7);
            p8 = populations(i,j,8);

            t = 1.0 / tau;  // TIME STEP IS IN LATTICE UNITS PORCO ZIO
            o_m_t = 1.0 - t;

            after_collision_populations(i,j,0) = o_m_t * p0 + t * equilibrium_populations(i,j,0);
            after_collision_populations(i,j,1) = o_m_t * p1 + t * equilibrium_populations(i,j,1);
            after_collision_populations(i,j,2) = o_m_t * p2 + t * equilibrium_populations(i,j,2);
            after_collision_populations(i,j,3) = o_m_t * p3 + t * equilibrium_populations(i,j,3);
            after_collision_populations(i,j,4) = o_m_t * p4 + t * equilibrium_populations(i,j,4);
            after_collision_populations(i,j,5) = o_m_t * p5 + t * equilibrium_populations(i,j,5);
            after_collision_populations(i,j,6) = o_m_t * p6 + t * equilibrium_populations(i,j,6);
            after_collision_populations(i,j,7) = o_m_t * p7 + t * equilibrium_populations(i,j,7);
            after_collision_populations(i,j,8) = o_m_t * p8 + t * equilibrium_populations(i,j,8);
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
        double p0, p1, p2, p3, p4, p5, p6, p7, p8;
        double t;
        double o_m_t;

        auto n_rows = populations.dimensions()[0];
        auto n_cols = populations.dimensions()[1];
        
        #pragma omp parallel for
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
     * @brief Most generic version of the BGK collision operator.
     * 
     * @tparam dim the total number of spatial dimensions.
     */
    template<std::size_t dim>
    class STDExecBGKCollisionPolicy : public CollisionPolicyTag, public STDExecTag
    {
    private:
        /// @brief Logger of the collision policy
        Logger l;

    public:
        STDExecBGKCollisionPolicy() 
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
    class STDExecBGKCollisionPolicy<2> : public CollisionPolicyTag, public STDExecTag
    {

    public:
        static double tau;
        
        /**
         * @brief Initializes the BGK collision policy with the
         * relaxation constant.
         * 
         * @param tau_ the relaxation constant
         */
        static void initialize(const double& tau_, const double& max_speed = 0.2)
        {
            tau = tau_;
            std::cout << "Inizialization done:" << std::endl;
            std::cout << " Tau: " << tau << std::endl;
            std::cout << " Max Speed: " << max_speed << std::endl;
            std::cout << "====================================================" << std::endl;
            if (tau <= 0.55)
            {
                std::cout << " The simulation might have stability issues." << std::endl;
                std::cout << " Theoretical LIMIT of tau:" << "tau > " << 0.5 + (1/8.0) * max_speed << std::endl;
            }
            std::cout << " Allowed Reynolds Number: " << max_speed * 3.0 / (tau - 0.5) << std::endl;
        }

        /**
         * @brief Inializes the value of the relaxation parameter based on the desired maximum velocity and Reynolds number
         * of the simulation. 
         * @note Quantities are expressed in lattice units.
         * 
         * @param max_velocity maximum velocity of the simulation 
         * @param reynolds Reynolds number of the simulation
         */
        static void allowed_tau(const double& max_speed, const double& reynolds, const double& speed_of_sound_squared = 1.0/3.0)
        {
            tau = 0.5 + max_speed / (speed_of_sound_squared * reynolds);
            std::cout << "Inizialization done:" << std::endl;
            std::cout << " Max Speed: " << max_speed << std::endl;
            std::cout << " Reynolds: " << reynolds << std::endl;
            std::cout << " Speed of Sound Squared: " << speed_of_sound_squared << std::endl;
            std::cout << "====================================================" << std::endl;
            std::cout << " Smallest theoretical Tau: " << tau << std::endl;
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
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double t;
            double o_m_t;

            //for(const auto& fluid_node : fluid_nodes)
            std::for_each(std::execution::par,fluid_nodes.begin(),fluid_nodes.end(),[&](const auto fluid_node)
            {
                i = fluid_node.coords[0];
                j = fluid_node.coords[1];

                p0 = populations(i,j,0);
                p1 = populations(i,j,1);
                p2 = populations(i,j,2);
                p3 = populations(i,j,3);
                p4 = populations(i,j,4);
                p5 = populations(i,j,5);
                p6 = populations(i,j,6);
                p7 = populations(i,j,7);
                p8 = populations(i,j,8);

                t = 1.0 / tau;  // TIME STEP IS IN LATTICE UNITS PORCO ZIO
                o_m_t = 1.0 - t;

                after_collision_populations(i,j,0) = o_m_t * p0 + t * equilibrium_populations(i,j,0);
                after_collision_populations(i,j,1) = o_m_t * p1 + t * equilibrium_populations(i,j,1);
                after_collision_populations(i,j,2) = o_m_t * p2 + t * equilibrium_populations(i,j,2);
                after_collision_populations(i,j,3) = o_m_t * p3 + t * equilibrium_populations(i,j,3);
                after_collision_populations(i,j,4) = o_m_t * p4 + t * equilibrium_populations(i,j,4);
                after_collision_populations(i,j,5) = o_m_t * p5 + t * equilibrium_populations(i,j,5);
                after_collision_populations(i,j,6) = o_m_t * p6 + t * equilibrium_populations(i,j,6);
                after_collision_populations(i,j,7) = o_m_t * p7 + t * equilibrium_populations(i,j,7);
                after_collision_populations(i,j,8) = o_m_t * p8 + t * equilibrium_populations(i,j,8);
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
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double t;
            double o_m_t;

            auto n_rows = populations.dimensions()[0];
            auto n_cols = populations.dimensions()[1];
            std::for_each(std::execution::par,fluid_nodes.begin(),fluid_nodes.end(),[&](const auto fluid_node)
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
        double p0, p1, p2, p3, p4, p5, p6, p7, p8;
        double t;
        double o_m_t;

        //for(const auto& open_boundary : open_boundary_nodes)
        std::for_each(std::execution::par,open_boundary_nodes.begin(),open_boundary_nodes.end(),[&](const auto open_boundary)
        {
            i = open_boundary.coords[0];
            j = open_boundary.coords[1];

            p0 = populations(i,j,0);
            p1 = populations(i,j,1);
            p2 = populations(i,j,2);
            p3 = populations(i,j,3);
            p4 = populations(i,j,4);
            p5 = populations(i,j,5);
            p6 = populations(i,j,6);
            p7 = populations(i,j,7);
            p8 = populations(i,j,8);

            t = 1.0 / tau;  // TIME STEP IS IN LATTICE UNITS PORCO ZIO
            o_m_t = 1.0 - t;

            after_collision_populations(i,j,0) = o_m_t * p0 + t * equilibrium_populations(i,j,0);
            after_collision_populations(i,j,1) = o_m_t * p1 + t * equilibrium_populations(i,j,1);
            after_collision_populations(i,j,2) = o_m_t * p2 + t * equilibrium_populations(i,j,2);
            after_collision_populations(i,j,3) = o_m_t * p3 + t * equilibrium_populations(i,j,3);
            after_collision_populations(i,j,4) = o_m_t * p4 + t * equilibrium_populations(i,j,4);
            after_collision_populations(i,j,5) = o_m_t * p5 + t * equilibrium_populations(i,j,5);
            after_collision_populations(i,j,6) = o_m_t * p6 + t * equilibrium_populations(i,j,6);
            after_collision_populations(i,j,7) = o_m_t * p7 + t * equilibrium_populations(i,j,7);
            after_collision_populations(i,j,8) = o_m_t * p8 + t * equilibrium_populations(i,j,8);
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
        double p0, p1, p2, p3, p4, p5, p6, p7, p8;
        double t;
        double o_m_t;

        auto n_rows = populations.dimensions()[0];
        auto n_cols = populations.dimensions()[1];
        std::for_each(std::execution::par,open_boundary_nodes.begin(),open_boundary_nodes.end(),[&](const auto open_boundary)
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
     * @brief Most generic version of the BGK collision operator.
     * 
     * @tparam dim the total number of spatial dimensions.
     */
    template<std::size_t dim>
    class OpenACCBGKCollisionPolicy : public CollisionPolicyTag, public OpenACCTag
    {
    private:
        /// @brief Logger of the collision policy
        Logger l;

    public:
        OpenACCBGKCollisionPolicy() 
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
    class OpenACCBGKCollisionPolicy<2> : public CollisionPolicyTag, public OpenACCTag
    {

    public:
        static double tau;
        
        /**
         * @brief Initializes the BGK collision policy with the
         * relaxation constant.
         * 
         * @param tau_ the relaxation constant
         */
        static void initialize(const double& tau_, const double& max_speed = 0.2)
        {
            tau = tau_;
            std::cout << "Inizialization done:" << std::endl;
            std::cout << " Tau: " << tau << std::endl;
            std::cout << " Max Speed: " << max_speed << std::endl;
            std::cout << "====================================================" << std::endl;
            if (tau <= 0.55)
            {
                std::cout << " The simulation might have stability issues." << std::endl;
                std::cout << " Theoretical LIMIT of tau:" << "tau > " << 0.5 + (1/8.0) * max_speed << std::endl;
            }
            std::cout << " Allowed Reynolds Number: " << max_speed * 3.0 / (tau - 0.5) << std::endl;
        }

        /**
         * @brief Inializes the value of the relaxation parameter based on the desired maximum velocity and Reynolds number
         * of the simulation. 
         * @note Quantities are expressed in lattice units.
         * 
         * @param max_velocity maximum velocity of the simulation 
         * @param reynolds Reynolds number of the simulation
         */
        static void allowed_tau(const double& max_speed, const double& reynolds, const double& speed_of_sound_squared = 1.0/3.0)
        {
            tau = 0.5 + max_speed / (speed_of_sound_squared * reynolds);
            std::cout << "Inizialization done:" << std::endl;
            std::cout << " Max Speed: " << max_speed << std::endl;
            std::cout << " Reynolds: " << reynolds << std::endl;
            std::cout << " Speed of Sound Squared: " << speed_of_sound_squared << std::endl;
            std::cout << "====================================================" << std::endl;
            std::cout << " Smallest theoretical Tau: " << tau << std::endl;
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
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double t;
            double o_m_t;
            
            #pragma acc parallel loop
            //for(const auto& fluid_node : fluid_nodes)
            for(size_t fnode=0; fnode < fluid_nodes.size(); fnode++)
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

                t = 1.0 / tau;  // TIME STEP IS IN LATTICE UNITS PORCO ZIO
                o_m_t = 1.0 - t;

                after_collision_populations(i,j,0) = o_m_t * p0 + t * equilibrium_populations(i,j,0);
                after_collision_populations(i,j,1) = o_m_t * p1 + t * equilibrium_populations(i,j,1);
                after_collision_populations(i,j,2) = o_m_t * p2 + t * equilibrium_populations(i,j,2);
                after_collision_populations(i,j,3) = o_m_t * p3 + t * equilibrium_populations(i,j,3);
                after_collision_populations(i,j,4) = o_m_t * p4 + t * equilibrium_populations(i,j,4);
                after_collision_populations(i,j,5) = o_m_t * p5 + t * equilibrium_populations(i,j,5);
                after_collision_populations(i,j,6) = o_m_t * p6 + t * equilibrium_populations(i,j,6);
                after_collision_populations(i,j,7) = o_m_t * p7 + t * equilibrium_populations(i,j,7);
                after_collision_populations(i,j,8) = o_m_t * p8 + t * equilibrium_populations(i,j,8);
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
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double t;
            double o_m_t;

            auto n_rows = populations.dimensions()[0];
            auto n_cols = populations.dimensions()[1];

            #pragma acc parallel loop
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
        double p0, p1, p2, p3, p4, p5, p6, p7, p8;
        double t;
        double o_m_t;

        #pragma acc parallel loop
        //for(const auto& open_boundary : open_boundary_nodes)
        for(size_t onode=0; onode < open_boundary_nodes.size(); onode++)
        {
            i = open_boundary_nodes[onode].coords[0];
            j = open_boundary_nodes[onode].coords[1];

            p0 = populations(i,j,0);
            p1 = populations(i,j,1);
            p2 = populations(i,j,2);
            p3 = populations(i,j,3);
            p4 = populations(i,j,4);
            p5 = populations(i,j,5);
            p6 = populations(i,j,6);
            p7 = populations(i,j,7);
            p8 = populations(i,j,8);

            t = 1.0 / tau;  // TIME STEP IS IN LATTICE UNITS PORCO ZIO
            o_m_t = 1.0 - t;

            after_collision_populations(i,j,0) = o_m_t * p0 + t * equilibrium_populations(i,j,0);
            after_collision_populations(i,j,1) = o_m_t * p1 + t * equilibrium_populations(i,j,1);
            after_collision_populations(i,j,2) = o_m_t * p2 + t * equilibrium_populations(i,j,2);
            after_collision_populations(i,j,3) = o_m_t * p3 + t * equilibrium_populations(i,j,3);
            after_collision_populations(i,j,4) = o_m_t * p4 + t * equilibrium_populations(i,j,4);
            after_collision_populations(i,j,5) = o_m_t * p5 + t * equilibrium_populations(i,j,5);
            after_collision_populations(i,j,6) = o_m_t * p6 + t * equilibrium_populations(i,j,6);
            after_collision_populations(i,j,7) = o_m_t * p7 + t * equilibrium_populations(i,j,7);
            after_collision_populations(i,j,8) = o_m_t * p8 + t * equilibrium_populations(i,j,8);
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
        double p0, p1, p2, p3, p4, p5, p6, p7, p8;
        double t;
        double o_m_t;

        auto n_rows = populations.dimensions()[0];
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
        }
    }

    };



    // initialization of the relaxation constant in the 2-D BGK collision operator.
    double BGKCollisionPolicy<2>::tau = 0.0;
    double OMPBGKCollisionPolicy<2>::tau = 0.0;
    double STDExecBGKCollisionPolicy<2>::tau = 0.0;
    double OpenACCBGKCollisionPolicy<2>::tau = 0.0;

}; // namespace llalbm::core::collisions

#endif // LLALBM_BGK_COLLISION_HPP