/**
 * @file PSBounceBackPolicy.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-02-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef LLALBM_PSBOUNCEBACKPOLICY_HPP
#define LLALBM_PSBOUNCEBACKPOLICY_HPP

#include <vector>
#include <map>
#include <array>
#include <execution>
// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "../../utils/loggers/Logger.hpp"
#include "../../utils/aliases.hpp" 
#include "../PolicyTypes.hpp"
// =======================================


namespace llalbm::core::boundaries
{
    
    using namespace llalbm::util::logger;
    using namespace llalbm::core;  
    //using namespace Eigen;  

    template<
        std::size_t dim>
    class PSBounceBackPolicy : public BoundaryPolicyTag, public SequentialTag{
    private:
        /* data */
        Logger log;
    public:
        PSBounceBackPolicy(/* args */)
        : log("Partially Saturated BB" , std::cout)
        {
            log.error("Partially Saturated BB not implemented");
            assert(false);
        }
    };

    template<>
    class PSBounceBackPolicy<2> : public BoundaryPolicyTag, public SequentialTag{
            
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

        static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
        {

            Eigen::Index i, j;
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;

            auto n_rows = populations.dimensions()[0];
            auto n_cols = populations.dimensions()[1];
            for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++) // * per castare il pointer????
            {
                i = boundary_coord[bnode].coords[0];
                j = boundary_coord[bnode].coords[1];

                p0 = populations(i, j, 0);
                p1 = populations(i, j, 1);
                p2 = populations(i, j, 2);
                p3 = populations(i, j, 3);
                p4 = populations(i, j, 4);
                p5 = populations(i, j, 5);
                p6 = populations(i, j, 6);
                p7 = populations(i, j, 7);
                p8 = populations(i, j, 8);

                // populations(i,j,0) = p0;

                // Compute the bounce back populations
                /*
                            j
                            6 2 5
                        i  3 0 1
                            7 4 8
                */
                
                // If the node is not on the left wall we can propagate to the left
                if(j != 0)
                {
                    populations(i, j-1, 3) = p1;
                }
                // If the node is not on the bottom wall we can propagate downwards
                if(i != n_rows-1)
                {
                    populations(i+1, j, 4) = p2;
                }
                // If the node is not on the right wall we can propagate to the right
                if(j != n_cols-1)
                {
                    populations(i, j+1, 1) = p3;
                }
                // If the node is not on the top wall we can propagate upwards
                if(i != 0)
                {
                    populations(i-1, j, 2) = p4;
                }
                // If the node is not on the bottom left corner we can propagate to the bottom left
                if(i != n_rows-1 && j != 0)
                {
                    populations(i+1, j-1, 7) = p5;
                }
                // If the node is not on the bottom right corner we can propagate to the bottom right
                if(i != n_rows-1 && j != n_cols-1)
                {
                    populations(i+1, j+1, 8) = p6;
                }
                // If the node is not on the top right corner we can propagate to the top right
                if(i != 0 && j != n_cols-1)
                {
                    populations(i-1, j+1, 5) = p7;
                }
                // If the node is not on the top left corner we can propagate to the top left
                if(i != 0 && j != 0)
                {
                    populations(i-1, j-1, 6) = p8;
                }
            }
        }

        static void update_boundaries(Tensor<double, 3> &populations, std::vector<ObstaclePoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
        {
            double rho,ux,uy,u2;
            constexpr double two9 = (2.0/9.0);
            constexpr double one18 = (1.0/18.0);
            constexpr double one36 = (1.0/36.0);
            double e0, e1, e2, e3, e4, e5, e6, e7, e8;
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double t;
            const std::size_t obstacles = boundary_coord.size();

            for (std::size_t o_node = 0; o_node < obstacles; ++o_node)
            {
                const Eigen::Index i = boundary_coord[o_node].coords[0];
                const Eigen::Index j = boundary_coord[o_node].coords[1];

                const double b = boundary_coord[o_node].get_b();

                // The first thing to do is define useful constant for the BGK method
                p1 = populations(i,j,1);
                p2 = populations(i,j,2);
                p3 = populations(i,j,3);
                p4 = populations(i,j,4);
                p5 = populations(i,j,5);
                p6 = populations(i,j,6);
                p7 = populations(i,j,7);
                p8 = populations(i,j,8);

                t = 1.0 / tau;

                const std::bitset<9> bb_mask = boundary_coord[o_node].directions;

                // Now, we calculate the effect of the boundary directly on the fluid node
                // For each useful fluid node, we take its density, its velocities and calculate the useful equilibriums
                if (bb_mask.test(1))
                {
                    rho = global_rho(i, j+1);
                    ux = global_u(i, j+1, 0);
                    uy = global_u(i, j+1, 1);
                    u2 = ux * ux + uy * uy;
                    e1 = one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                    e3 = one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                    populations(i, j+1, 1) = p3 - (1.0 - b) * (p3 - e3) * t + b * ((populations(i, j+1, 1) - e1) - (p3 - 2 * one18 * global_rho(i,j)));
                }
                if (bb_mask.test(2))
                {
                    rho = global_rho(i-1,j);
                    ux = global_u(i-1,j,0);
                    uy = global_u(i-1,j,1);
                    u2 = ux * ux + uy * uy;
                    e2 = one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                    e4 = one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                    populations(i-1, j, 2) = p4 - (1.0 - b) * (p4 - e4) * t + b * ((populations(i-1, j, 2) - e2) - (p4 - 2 * one18 * global_rho(i,j)));
                }
                if (bb_mask.test(3))
                {
                    rho = global_rho(i,j-1);
                    ux = global_u(i,j-1,0);
                    uy = global_u(i,j-1,1);
                    u2 = ux * ux + uy * uy;
                    e1 = one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                    e3 = one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                    populations(i, j-1, 3) = p1 - (1.0 - b) * (p1 - e1) * t + b * ((populations(i, j-1, 3) - e3) - (p1 - 2 * one18 * global_rho(i,j)));
                }
                if (bb_mask.test(4))
                {
                    rho = global_rho(i+1,j);
                    ux = global_u(i+1,j,0);
                    uy = global_u(i+1,j,1);
                    u2 = ux * ux + uy * uy;
                    e2 = one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                    e4 = one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                    populations(i+1, j, 4) = p2 - (1.0 - b) * (p2 - e2) * t + b * ((populations(i+1, j, 4) - e4) - (p2 - 2 * one18 * global_rho(i,j)));
                }
                if (bb_mask.test(5))
                {
                    rho = global_rho(i-1,j+1);
                    ux = global_u(i-1,j+1,0);
                    uy = global_u(i-1,j+1,1);
                    u2 = ux * ux + uy * uy;
                    e5 = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    e7 = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    populations(i-1, j+1, 5) = p7 - (1.0 - b) * (p7 - e7) * t + b * ((populations(i-1, j+1, 5) - e5) - (p7 - one36 * global_rho(i,j)));
                }
                if (bb_mask.test(6))
                {
                    rho = global_rho(i-1,j-1);
                    ux = global_u(i-1,j-1,0);
                    uy = global_u(i-1,j-1,1);
                    u2 = ux * ux + uy * uy;
                    e6 = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    e8 = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    populations(i-1, j-1, 6) = p8 - (1.0 - b) * (p8 - e8) * t + b * ((populations(i-1, j-1, 6) - e6) - (p8 - one36 * global_rho(i,j)));
                }
                if (bb_mask.test(7))
                {
                    rho = global_rho(i+1,j-1);
                    ux = global_u(i+1,j-1,0);
                    uy = global_u(i+1,j-1,1);
                    u2 = ux * ux + uy * uy;
                    e5 = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    e7 = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    populations(i+1, j-1, 7) = p5 - (1.0 - b) * (p5 - e5) * t + b * ((populations(i+1, j-1, 7) - e7) - (p5 - one36 * global_rho(i,j)));
                }
                if (bb_mask.test(8))
                {
                    rho = global_rho(i+1,j+1);
                    ux = global_u(i+1,j+1,0);
                    uy = global_u(i+1,j+1,1);
                    u2 = ux * ux + uy * uy;
                    e6 = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    e8 = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    populations(i+1, j+1, 8) = p6 - (1.0 - b) * (p6 - e6) * t + b * ((populations(i+1, j+1, 8) - e8) - (p6 - one36 * global_rho(i,j)));
                }
                
            }
        }
    };


    template<
        std::size_t dim>
    class OMPPSBounceBackPolicy : public BoundaryPolicyTag, public OMPTag{
    private:
        /* data */
        Logger log;
    public:
        OMPPSBounceBackPolicy(/* args */)
        : log("Partially Saturated BB" , std::cout)
        {
            log.error("Partially Saturated BB not implemented");
            assert(false);
        }
    };

    template<>
    class OMPPSBounceBackPolicy<2> : public BoundaryPolicyTag, public OMPTag{
            
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

        static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
        {

            Eigen::Index i, j;
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;

            auto n_rows = populations.dimensions()[0];
            auto n_cols = populations.dimensions()[1];
            #pragma omp parallel for private(i,j,p0,p1,p2,p3,p4,p5,p6,p7,p8, n_rows, n_cols)
            for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++) // * per castare il pointer????
            {
                i = boundary_coord[bnode].coords[0];
                j = boundary_coord[bnode].coords[1];

                p0 = populations(i, j, 0);
                p1 = populations(i, j, 1);
                p2 = populations(i, j, 2);
                p3 = populations(i, j, 3);
                p4 = populations(i, j, 4);
                p5 = populations(i, j, 5);
                p6 = populations(i, j, 6);
                p7 = populations(i, j, 7);
                p8 = populations(i, j, 8);

                // populations(i,j,0) = p0;

                // Compute the bounce back populations
                /*
                            j
                            6 2 5
                        i  3 0 1
                            7 4 8
                */
                
                // If the node is not on the left wall we can propagate to the left
                if(j != 0)
                {
                    populations(i, j-1, 3) = p1;
                }
                // If the node is not on the bottom wall we can propagate downwards
                if(i != n_rows-1)
                {
                    populations(i+1, j, 4) = p2;
                }
                // If the node is not on the right wall we can propagate to the right
                if(j != n_cols-1)
                {
                    populations(i, j+1, 1) = p3;
                }
                // If the node is not on the top wall we can propagate upwards
                if(i != 0)
                {
                    populations(i-1, j, 2) = p4;
                }
                // If the node is not on the bottom left corner we can propagate to the bottom left
                if(i != n_rows-1 && j != 0)
                {
                    populations(i+1, j-1, 7) = p5;
                }
                // If the node is not on the bottom right corner we can propagate to the bottom right
                if(i != n_rows-1 && j != n_cols-1)
                {
                    populations(i+1, j+1, 8) = p6;
                }
                // If the node is not on the top right corner we can propagate to the top right
                if(i != 0 && j != n_cols-1)
                {
                    populations(i-1, j+1, 5) = p7;
                }
                // If the node is not on the top left corner we can propagate to the top left
                if(i != 0 && j != 0)
                {
                    populations(i-1, j-1, 6) = p8;
                }
            }
        }

        static void update_boundaries(Tensor<double, 3> &populations, std::vector<ObstaclePoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
        {
            double rho,ux,uy,u2;
            constexpr double two9 = (2.0/9.0);
            constexpr double one18 = (1.0/18.0);
            constexpr double one36 = (1.0/36.0);
            double e0, e1, e2, e3, e4, e5, e6, e7, e8;
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double t;
            const std::size_t obstacles = boundary_coord.size();

            #pragma omp parallel for private(p0,p1,p2,p3,p4,p5,p6,p7,p8,e0,e1,e2,e3,e4,e5,e6,e7,e8,rho,ux,uy,u2,t)
            for (std::size_t o_node = 0; o_node < obstacles; ++o_node)
            {
                const Eigen::Index i = boundary_coord[o_node].coords[0];
                const Eigen::Index j = boundary_coord[o_node].coords[1];

                const double b = boundary_coord[o_node].get_b();

                // The first thing to do is define useful constant for the BGK method
                p1 = populations(i,j,1);
                p2 = populations(i,j,2);
                p3 = populations(i,j,3);
                p4 = populations(i,j,4);
                p5 = populations(i,j,5);
                p6 = populations(i,j,6);
                p7 = populations(i,j,7);
                p8 = populations(i,j,8);

                t = 1.0 / tau;

                const std::bitset<9> bb_mask = boundary_coord[o_node].directions;

                // Now, we calculate the effect of the boundary directly on the fluid node
                // For each useful fluid node, we take its density, its velocities and calculate the useful equilibriums
                if (bb_mask.test(1))
                {
                    rho = global_rho(i, j+1);
                    ux = global_u(i, j+1, 0);
                    uy = global_u(i, j+1, 1);
                    u2 = ux * ux + uy * uy;
                    e1 = one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                    e3 = one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                    populations(i, j+1, 1) = p3 - (1.0 - b) * (p3 - e3) * t + b * ((populations(i, j+1, 1) - e1) - (p3 - 2 * one18 * global_rho(i,j)));
                }
                if (bb_mask.test(2))
                {
                    rho = global_rho(i-1,j);
                    ux = global_u(i-1,j,0);
                    uy = global_u(i-1,j,1);
                    u2 = ux * ux + uy * uy;
                    e2 = one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                    e4 = one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                    populations(i-1, j, 2) = p4 - (1.0 - b) * (p4 - e4) * t + b * ((populations(i-1, j, 2) - e2) - (p4 - 2 * one18 * global_rho(i,j)));
                }
                if (bb_mask.test(3))
                {
                    rho = global_rho(i,j-1);
                    ux = global_u(i,j-1,0);
                    uy = global_u(i,j-1,1);
                    u2 = ux * ux + uy * uy;
                    e1 = one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                    e3 = one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                    populations(i, j-1, 3) = p1 - (1.0 - b) * (p1 - e1) * t + b * ((populations(i, j-1, 3) - e3) - (p1 - 2 * one18 * global_rho(i,j)));
                }
                if (bb_mask.test(4))
                {
                    rho = global_rho(i+1,j);
                    ux = global_u(i+1,j,0);
                    uy = global_u(i+1,j,1);
                    u2 = ux * ux + uy * uy;
                    e2 = one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                    e4 = one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                    populations(i+1, j, 4) = p2 - (1.0 - b) * (p2 - e2) * t + b * ((populations(i+1, j, 4) - e4) - (p2 - 2 * one18 * global_rho(i,j)));
                }
                if (bb_mask.test(5))
                {
                    rho = global_rho(i-1,j+1);
                    ux = global_u(i-1,j+1,0);
                    uy = global_u(i-1,j+1,1);
                    u2 = ux * ux + uy * uy;
                    e5 = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    e7 = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    populations(i-1, j+1, 5) = p7 - (1.0 - b) * (p7 - e7) * t + b * ((populations(i-1, j+1, 5) - e5) - (p7 - one36 * global_rho(i,j)));
                }
                if (bb_mask.test(6))
                {
                    rho = global_rho(i-1,j-1);
                    ux = global_u(i-1,j-1,0);
                    uy = global_u(i-1,j-1,1);
                    u2 = ux * ux + uy * uy;
                    e6 = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    e8 = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    populations(i-1, j-1, 6) = p8 - (1.0 - b) * (p8 - e8) * t + b * ((populations(i-1, j-1, 6) - e6) - (p8 - one36 * global_rho(i,j)));
                }
                if (bb_mask.test(7))
                {
                    rho = global_rho(i+1,j-1);
                    ux = global_u(i+1,j-1,0);
                    uy = global_u(i+1,j-1,1);
                    u2 = ux * ux + uy * uy;
                    e5 = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    e7 = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    populations(i+1, j-1, 7) = p5 - (1.0 - b) * (p5 - e5) * t + b * ((populations(i+1, j-1, 7) - e7) - (p5 - one36 * global_rho(i,j)));
                }
                if (bb_mask.test(8))
                {
                    rho = global_rho(i+1,j+1);
                    ux = global_u(i+1,j+1,0);
                    uy = global_u(i+1,j+1,1);
                    u2 = ux * ux + uy * uy;
                    e6 = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    e8 = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    populations(i+1, j+1, 8) = p6 - (1.0 - b) * (p6 - e6) * t + b * ((populations(i+1, j+1, 8) - e8) - (p6 - one36 * global_rho(i,j)));
                }
                
            }
        }
    };

    template<
        std::size_t dim>
    class STDExecPSBounceBackPolicy : public BoundaryPolicyTag, public STDExecTag{
    private:
        /* data */
        Logger log;
    public:
        STDExecPSBounceBackPolicy(/* args */)
        : log("Partially Saturated BB" , std::cout)
        {
            log.error("Partially Saturated BB not implemented");
            assert(false);
        }
    };

    template<>
    class STDExecPSBounceBackPolicy<2> : public BoundaryPolicyTag, public STDExecTag{
            
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

        static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
        {

            Eigen::Index i, j;
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;

            auto n_rows = populations.dimensions()[0];
            auto n_cols = populations.dimensions()[1];
            //for (size_t bnode = 0; bnode < boundary_coord.size(); bnode++) // * per castare il pointer????
            std::for_each(std::execution::par, boundary_coord.begin(), boundary_coord.end(), [&](const auto bnode)
            {
                i = bnode.coords[0];
                j = bnode.coords[1];

                p0 = populations(i, j, 0);
                p1 = populations(i, j, 1);
                p2 = populations(i, j, 2);
                p3 = populations(i, j, 3);
                p4 = populations(i, j, 4);
                p5 = populations(i, j, 5);
                p6 = populations(i, j, 6);
                p7 = populations(i, j, 7);
                p8 = populations(i, j, 8);

                // populations(i,j,0) = p0;

                // Compute the bounce back populations
                /*
                            j
                            6 2 5
                        i  3 0 1
                            7 4 8
                */
                
                // If the node is not on the left wall we can propagate to the left
                if(j != 0)
                {
                    populations(i, j-1, 3) = p1;
                }
                // If the node is not on the bottom wall we can propagate downwards
                if(i != n_rows-1)
                {
                    populations(i+1, j, 4) = p2;
                }
                // If the node is not on the right wall we can propagate to the right
                if(j != n_cols-1)
                {
                    populations(i, j+1, 1) = p3;
                }
                // If the node is not on the top wall we can propagate upwards
                if(i != 0)
                {
                    populations(i-1, j, 2) = p4;
                }
                // If the node is not on the bottom left corner we can propagate to the bottom left
                if(i != n_rows-1 && j != 0)
                {
                    populations(i+1, j-1, 7) = p5;
                }
                // If the node is not on the bottom right corner we can propagate to the bottom right
                if(i != n_rows-1 && j != n_cols-1)
                {
                    populations(i+1, j+1, 8) = p6;
                }
                // If the node is not on the top right corner we can propagate to the top right
                if(i != 0 && j != n_cols-1)
                {
                    populations(i-1, j+1, 5) = p7;
                }
                // If the node is not on the top left corner we can propagate to the top left
                if(i != 0 && j != 0)
                {
                    populations(i-1, j-1, 6) = p8;
                }
            });
        }

        static void update_boundaries(Tensor<double, 3> &populations, std::vector<ObstaclePoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
        {
            double rho,ux,uy,u2;
            constexpr double two9 = (2.0/9.0);
            constexpr double one18 = (1.0/18.0);
            constexpr double one36 = (1.0/36.0);
            double e0, e1, e2, e3, e4, e5, e6, e7, e8;
            double p0, p1, p2, p3, p4, p5, p6, p7, p8;
            double t;
            //const std::size_t obstacles = boundary_coord.size();

            //for (std::size_t o_node = 0; o_node < obstacles; ++o_node)
            std::for_each(std::execution::par, boundary_coord.begin(), boundary_coord.end(), [&](const auto o_node)
            {
                const Eigen::Index i = o_node.coords[0];
                const Eigen::Index j = o_node.coords[1];

                const double b = o_node.get_b();

                // The first thing to do is define useful constant for the BGK method
                p1 = populations(i,j,1);
                p2 = populations(i,j,2);
                p3 = populations(i,j,3);
                p4 = populations(i,j,4);
                p5 = populations(i,j,5);
                p6 = populations(i,j,6);
                p7 = populations(i,j,7);
                p8 = populations(i,j,8);

                t = 1.0 / tau;

                const std::bitset<9> bb_mask = o_node.directions;

                // Now, we calculate the effect of the boundary directly on the fluid node
                // For each useful fluid node, we take its density, its velocities and calculate the useful equilibriums
                if (bb_mask.test(1))
                {
                    rho = global_rho(i, j+1);
                    ux = global_u(i, j+1, 0);
                    uy = global_u(i, j+1, 1);
                    u2 = ux * ux + uy * uy;
                    e1 = one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                    e3 = one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                    populations(i, j+1, 1) = p3 - (1.0 - b) * (p3 - e3) * t + b * ((populations(i, j+1, 1) - e1) - (p3 - 2 * one18 * global_rho(i,j)));
                }
                if (bb_mask.test(2))
                {
                    rho = global_rho(i-1,j);
                    ux = global_u(i-1,j,0);
                    uy = global_u(i-1,j,1);
                    u2 = ux * ux + uy * uy;
                    e2 = one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                    e4 = one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                    populations(i-1, j, 2) = p4 - (1.0 - b) * (p4 - e4) * t + b * ((populations(i-1, j, 2) - e2) - (p4 - 2 * one18 * global_rho(i,j)));
                }
                if (bb_mask.test(3))
                {
                    rho = global_rho(i,j-1);
                    ux = global_u(i,j-1,0);
                    uy = global_u(i,j-1,1);
                    u2 = ux * ux + uy * uy;
                    e1 = one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                    e3 = one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                    populations(i, j-1, 3) = p1 - (1.0 - b) * (p1 - e1) * t + b * ((populations(i, j-1, 3) - e3) - (p1 - 2 * one18 * global_rho(i,j)));
                }
                if (bb_mask.test(4))
                {
                    rho = global_rho(i+1,j);
                    ux = global_u(i+1,j,0);
                    uy = global_u(i+1,j,1);
                    u2 = ux * ux + uy * uy;
                    e2 = one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                    e4 = one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                    populations(i+1, j, 4) = p2 - (1.0 - b) * (p2 - e2) * t + b * ((populations(i+1, j, 4) - e4) - (p2 - 2 * one18 * global_rho(i,j)));
                }
                if (bb_mask.test(5))
                {
                    rho = global_rho(i-1,j+1);
                    ux = global_u(i-1,j+1,0);
                    uy = global_u(i-1,j+1,1);
                    u2 = ux * ux + uy * uy;
                    e5 = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    e7 = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    populations(i-1, j+1, 5) = p7 - (1.0 - b) * (p7 - e7) * t + b * ((populations(i-1, j+1, 5) - e5) - (p7 - one36 * global_rho(i,j)));
                }
                if (bb_mask.test(6))
                {
                    rho = global_rho(i-1,j-1);
                    ux = global_u(i-1,j-1,0);
                    uy = global_u(i-1,j-1,1);
                    u2 = ux * ux + uy * uy;
                    e6 = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    e8 = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    populations(i-1, j-1, 6) = p8 - (1.0 - b) * (p8 - e8) * t + b * ((populations(i-1, j-1, 6) - e6) - (p8 - one36 * global_rho(i,j)));
                }
                if (bb_mask.test(7))
                {
                    rho = global_rho(i+1,j-1);
                    ux = global_u(i+1,j-1,0);
                    uy = global_u(i+1,j-1,1);
                    u2 = ux * ux + uy * uy;
                    e5 = one36 * rho * (1 + 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    e7 = one36 * rho * (1 - 3 * (ux + uy) + 9 * ux * uy + 3 * u2);
                    populations(i+1, j-1, 7) = p5 - (1.0 - b) * (p5 - e5) * t + b * ((populations(i+1, j-1, 7) - e7) - (p5 - one36 * global_rho(i,j)));
                }
                if (bb_mask.test(8))
                {
                    rho = global_rho(i+1,j+1);
                    ux = global_u(i+1,j+1,0);
                    uy = global_u(i+1,j+1,1);
                    u2 = ux * ux + uy * uy;
                    e6 = one36 * rho * (1 - 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    e8 = one36 * rho * (1 + 3 * (ux - uy) - 9 * ux * uy + 3 * u2);
                    populations(i+1, j+1, 8) = p6 - (1.0 - b) * (p6 - e6) * t + b * ((populations(i+1, j+1, 8) - e8) - (p6 - one36 * global_rho(i,j)));
                }
                
            });
        }
    };



    // initialization of the relaxation constant in the 2-D PS Bounce-Back.
    double PSBounceBackPolicy<2>::tau = 0.0;

    // initialization of the relaxation constant in the 2-D OMP PS Bounce-Back.
    double OMPPSBounceBackPolicy<2>::tau = 0.0;

    double STDExecPSBounceBackPolicy<2>::tau = 0.0;
}

#endif