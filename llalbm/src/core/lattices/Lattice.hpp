#ifndef LLALBM_LATTICE_HPP
#define LLALBM_LATTICE_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <map>
#include <array>
#include <cassert>
// ======================================

// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "Logger.hpp"
// =======================================

namespace llalbm::core
{
    using namespace Eigen;
    using namespace llalbm::util::logger;

    /// Generic d-dimensional point described by its set of coordinates
    template<std::size_t d>
    using Point = Matrix<std::size_t, d, 1>;

    // TODO: partial specialization for a distributed lattice
    /**
     * The lattice class represents one of the most important components of a Lattice-Bolzmann simulation: the computational domain
     * onto which the simulation is performed. The class makes extensive use of the Eigen Linear Algebra library: the grid of points as
     * well as the particle distributions are stored as Eigen tensors, and are built based on the supplied spatial information when the object
     * is first instantiated. Templating allows for great versatility: with just a few limitations (mainly when using a GPU) the user is able to 
     * decide which policies to use or is free to implement other policies for the specific use-case.
     * 
    */
    template<
        typename ParallelizationPolicy  ,   // Policy managing the parallelization of the lattice and of the algorithm
        std::size_t dim                 ,   // Spacial dimensions of the simulation.
        typename CollisionPolicy        ,   // Policy managing the interaction between two fluid nodes.
        typename BoundaryPolicy         ,   // Policy managing the interaction between a fluid node and a boundary node
        typename InitializationPolicy       // Policy managing the initialization of the lattice.
    >
    class Lattice
    {
    private:

        /// @brief Policy used to manage the interaction between two fluid nodes.
        CollisionPolicy collision_policy;

        /// @brief Policy used to manage the interaction between a fluid node and a boundary node.
        BoundaryPolicy boundary_policy;

        /// @brief Policy used to manage the initialization of the lattice nodes (both when the lattice is being created)
        /// and at the beginning of each iteration to set the values of boundary nodes.
        InitializationPolicy initialization_policy;

        /// @brief Population tensor: at coordinates x (with x being an n-dimensional vector), we store
        /// the total vector of populations (hence the dim + 1)
        Tensor<double, dim + 1> populations;

        /// @brief Equilibrium populations tensor.
        Tensor<double, dim + 1> equilibrium_populations;

        // TODO: undestrand whether it would be better to move this tensor directly inside the policy.
        /// @brief Tensor of populations after the collision step and before the streaming step.
        Tensor<double, dim + 1> after_collision_populations;

        /// @brief Global density tensor.
        Tensor<double, dim> global_rho;

        /// @brief Global velocity tensor.
        Tensor<double, dim + 1> global_u;

        /// @brief Number of elements per lattice dimension.
        const std::array<Eigen::Index, dim> lattice_dimensions; 

        /// @brief List of coordinates of boundary_nodes.
        std::vector<Point<dim>> boundary_nodes;

        /// @brief List of coordinates of the edges of obstacle nodes.
        std::vector<Point<dim>> obstacle_nodes;

        /// @brief Number of velocities in the velocity set
        const std::size_t q;

        /// @brief Logger.
        Logger logger;

    public:
        
        /**
         * @brief Constructs a new Lattice object by instantiating the data structures and by initializing
         * the policies with useful information.
         * 
         * @param lattice_dimensions_ array indicating the extension of each spatial dimension in terms of grid points.
         * @param q_ number of velocities.
         * @param out_stream output stream used for logging.
         */
        Lattice(
            const std::array<Eigen::Index, dim>& lattice_dimensions_,
            const std::size_t& q_,
            std::ostream& out_stream = std::cout
        )
        :   lattice_dimensions (lattice_dimensions_)
        ,   q (q)
        ,   Logger("LATTICE", out_stream)
        {
            // Checking validity of input data
            assert(lattice_dimensions.size() == dim && "[ERROR] The number of dimension data in the Lattice constructor"+
            " must be equal to the number of spatial dimensions.");

            logger.info("Constructing Lattice Object...");
            logger.info("Number of dimensions = " + dim);
            logger.info("Numero of velocities = " + q);
        
            // Building the population tensor
            std::array<Eigen::Index, dim+1> populations_data;
            for (std::size_t i = 0; i < dim; ++i) {
               populations_data[i] = lattice_dimensions[i];
            }
            populations_data[dim] = q;
            
            populations = Tensor<double, dim + 1>(populations_data);
            equilibrium_populations = Tensor<double, dim + 1>(populations_data);
            after_collision_populations = Tensor<double, dim + 1>(populations_data);
            logger.info("Populations tensors built");

            // Building the global density tensor
            global_rho = Tensor<double, dim>(lattice_dimensions);
            logger.info("Global density tensor built");

            // Building the gloabl velocity tensor
            std::array<Eigen::Index, dim + 1> global_u_data;
            for (std::size_t i = 0; i < dim; ++i) {
               global_u_data[i] = lattice_dims[i];
            }
            global_u_data[dim] = dim;
            global_u = Tensor<double, dim + 1>(global_u_data);
            logger.info("Global velocity tensor build");

            //TODO: once collisions/boundary etc have been defined, attach relevant data.
            //TODO:
            //TODO:

            // Then, lattice construction is complete.
            logger.info("Lattice is ready.")
        }

    };
} // namespace llalbm::core

#endif // LLALBM_LATTICE_HPP