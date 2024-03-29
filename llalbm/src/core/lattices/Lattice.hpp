#ifndef LLALBM_LATTICE_HPP
#define LLALBM_LATTICE_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
// ======================================

// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "../../utils/loggers/Logger.hpp"
#include "../../utils/aliases.hpp"
#include "../../utils/writers/LatticeWriter.hpp"
#include "../../utils/readers/LatticeReader.hpp"
#include "../equilibriums/Equilibrium.hpp"
// =======================================

//6 2 5
//3 0 1
//7 4 8

namespace llalbm::core
{
    //using namespace Eigen;
    using namespace llalbm::util::logger;

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
        std::size_t dim                 ,   // Spacial dimensions of the simulation.
        typename CollisionPolicy        ,   // Policy managing the interaction between two fluid nodes.
        typename WallPolicy             ,   // Policy managing the interaction between a fluid node and a wall node.
        typename ObstaclePolicy         ,   // Policy managing the interaction between a fluid node and an internal obstacle.
        typename InletPolicy            ,   // Policy managing the interaction between an internal fluid node and an inlet node.
        typename OutletPolicy           ,   // Policy managing the interaction between an internal fluid node and an outlet policy.s
        typename InitializationPolicy       // Policy managing the initialization of the lattice.
    >
    class Lattice
    {
    private:

        // ========= POLICIES =========

        /// @brief Policy used to manage the interaction between two fluid nodes.
        CollisionPolicy collision_policy;

        /// @brief Policy used to manage the interaction between a fluid node and a boundary node.
        WallPolicy boundary_policy;

        /// @brief Policy used to manage the interaction between a fluid node and an internal obstacle.
        ObstaclePolicy obstacle_policy;

        /// @brief Policy used to manage the interaction between a fluid node and an inlet node.
        InletPolicy inlet_policy;

        /// @brief Policy used to manage the interaction between a fluid node and an outlet node
        OutletPolicy outlet_policy;

        /// @brief Policy used to manage the initialization of the lattice nodes (both when the lattice is being created)
        /// and at the beginning of each iteration to set the values of boundary nodes.
        InitializationPolicy initialization_policy;

        // ========= TENSORS OF THE LATTICE =========

        /// @brief Population tensor: at coordinates x (with x being an n-dimensional vector), we store
        /// the total vector of populations (hence the dim + 1)
        Tensor<double, dim + 1> populations;

        /// @brief Equilibrium populations tensor.
        Tensor<double, dim + 1> equilibrium_populations;

        /// @brief Tensor of populations after the collision step and before the streaming step.
        Tensor<double, dim + 1> after_collision_populations;

        /// @brief Global density tensor.
        Tensor<double, dim> global_rho;

        /// @brief Global velocity tensor.
        Tensor<double, dim + 1> global_u;

        /// @brief Number of elements per lattice dimension.
        std::array<Eigen::Index, dim> lattice_dimensions; 

        // ========= COORDINATES OF ALL NODES =========
        
        /// @brief List of coordinates of fluid nodes
        std::vector<Point<dim>> fluid_nodes;

        /// @brief List of coordinates of boundary nodes
        std::vector<BoundaryPoint<dim>> boundary_coord;

        /// @brief List of coordinates of inlet nodes.
        std::vector<BoundaryPoint<dim>> inlet_nodes_coord;

        /// @brief List of coordinates of outlet nodes.
        std::vector<BoundaryPoint<dim>> outlet_nodes_coord;

        /// @brief List of coordinates of the edges of obstacle nodes.
        std::vector<BoundaryPoint<dim>> obstacle_nodes;

        /// @brief Number of velocities in the velocity set
        const std::size_t q;

        // ========= MISCELLANEOUS =========

        /// @brief Logger.
        Logger logger;

        /// @brief Boolean used to save macroscopic quantities
        bool save;
    
    private:

        // ========================================================================================= 
        //                        GATEKEEPING FUNCTIONS FOR LATTICE CONSTRUCTION
        // ========================================================================================= 

        /**
         * @brief Initializes the velocity sets for the D2Q9 and D3Q19 lattices.
         */
        void initialize_velocity_sets()
        {

            D2Q9 << 0 , 0, 4/9,
                1 , 0, 1/9,
                0 , 1, 1/9,
                -1, 0, 1/9,
                0 ,-1, 1/9,
                1 , 1, 1/36,
                -1, 1, 1/36,
                -1,-1, 1/36,
                1 ,-1, 1/36;
            
            D3Q19 << 0, 0, 0, 1/3,
                1, 0, 0, 1/18,
                -1, 0, 0, 1/18,
                0, 1, 0, 1/18,
                0, -1, 0, 1/18,
                0, 0, 1, 1/18,
                0, 0, -1, 1/18,
                1, 1, 0, 1/36,
                -1, -1, 0, 1/36,
                1, 0, 1, 1/36,
                -1, 0, -1, 1/36,
                0, 1, 1, 1/36,
                0, -1, -1, 1/36,
                1, -1, 0, 1/36,
                -1, 1, 0, 1/36,
                1, 0, -1, 1/36,
                -1, 0, 1, 1/36,
                0, 1, -1, 1/36,
                0, -1, 1, 1/36;
        }

        /**
         * @brief Builds the tensors that describe the lattice.
         */
        void build_tensors()
        {
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
               global_u_data[i] = lattice_dimensions[i];
            }
            global_u_data[dim] = dim;
            global_u = Tensor<double, dim + 1>(global_u_data);
            logger.info("Global velocity tensor build");
        }

        // =========================================================================================
        // =========================================================================================
        // =========================================================================================


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
        ,   q (q_)
        ,   logger("LATTICE", out_stream)
        {
            // Checking validity of input data
            assert(lattice_dimensions.size() == dim && 
            "[ERROR] The number of dimension data in the Lattice constructor \n must be equal to the number of spatial dimensions.");

            logger.info("Constructing Lattice Object...");
            logger.info("Number of dimensions = " + std::to_string(dim));
            logger.info("Number of velocities = " + std::to_string(q));


            // Build the tensors
            build_tensors();

            // Initializes the velocity sets
            initialize_velocity_sets();

            logger.info("Velocity Sets initialized");
            //TODO: once collisions/boundary etc have been defined, attach relevant data.
            //TODO: understand how to handle output directory creation (if here or in another class)

            // Attach the nodes vectors to the initialization policy
            initialization_policy.attach_nodes(
                inlet_nodes_coord,
                outlet_nodes_coord
            );
            logger.info("Nodes attached to the initialization policy");

            // Then, lattice construction is complete.
            logger.info("Lattice is ready.");
        }

        Lattice(
            const std::string& path_,
            const std::size_t& q_,
            std::ostream& out_stream = std::cout
        )
        :   q (q_)
        ,   logger("LATTICE", out_stream)
        {
            llalbm::util::reader::read_lattice_file(path_,fluid_nodes,boundary_coord,inlet_nodes_coord,outlet_nodes_coord,obstacle_nodes,lattice_dimensions);
            // Checking validity of input data
            assert(lattice_dimensions.size() == dim && 
            "[ERROR] The number of dimension data in the Lattice constructor \n must be equal to the number of spatial dimensions.");

            logger.info("Constructing Lattice Object...");
            logger.info("Number of dimensions = " + std::to_string(dim));
            logger.info("Number of velocities = " + std::to_string(q_));

            // Build the tensors
            build_tensors();

            // Initialize the velocity sets
            initialize_velocity_sets();
            

            logger.info("Velocity Sets initialized");

            // Attach the nodes vectors to the initialization policy
            initialization_policy.attach_nodes(
                inlet_nodes_coord,
                outlet_nodes_coord
            );
            logger.info("Nodes attached to the initialization policy");

            // Then, lattice construction is complete.
            logger.info("Lattice is ready.");
        }

        void perform_lbm(const double time, const double time_step = 0.05, const std::size_t save_step = 5, const bool should_save = true)
        {
            assert(time_step > 0 && "ERROR: time step must be greater than 0");
            assert(time > 0 && "ERROR: simulation time must be greater than 0");
            assert(time > time_step && "ERROR: simulation time must be greater than the time step");

            // Initialize the values of the tensors
            global_rho.setConstant(1);
            global_u.setConstant(0);

            // The output directory is created only if the user wants to save the results
            if (should_save)
            {
                llalbm::util::writer::create_output_directory();
            }
            else
            {
                logger.info("Results will not be saved");
            }
            
            // Compute the total number of iterations
            const std::size_t n_steps = time/time_step;
            std::size_t saved_file = 0;

            llalbm::core::equilibrium::Equilibrium<dim>::calc_equilibrium(fluid_nodes, populations, global_u, global_rho);
            for (std::size_t i = 0; i < n_steps; i++)
            {
                initialization_policy.update_nodes(i, global_u, global_rho);

                save = (i%save_step == 0 && should_save);

                // 1. The equilibrium populations are calculated for each node
                llalbm::core::equilibrium::Equilibrium<dim>::calc_equilibrium(fluid_nodes,equilibrium_populations,global_u,global_rho);
                
                //2. Perform the collisions
                collision_policy.collide(populations, equilibrium_populations, after_collision_populations, fluid_nodes, global_rho, global_u);

                //3. Streaming
                collision_policy.stream(populations, after_collision_populations, fluid_nodes);

                //4. Perform the collision at the boundaries
                boundary_policy.update_boundaries(populations, boundary_coord, global_rho, global_u);
                obstacle_policy.update_boundaries(populations, obstacle_nodes, global_rho, global_u);
                inlet_policy.update_boundaries(populations, inlet_nodes_coord, global_rho, global_u);
                outlet_policy.update_boundaries(populations, outlet_nodes_coord, global_rho, global_u);
                
                // 5. Save ux and uy on files
                if(save)
                {
                    llalbm::util::writer::write_lattice_file<dim>(global_u, saved_file);
                    saved_file++;
                }
            }
        }

        // ========================================================================================= 
        //                                     GETTERS AND SETTERS     
        // ========================================================================================= 

        /**
         * @brief Get the global u object
         * 
         * @return Eigen::Tensor<double, dim + 1>& 
         */
        Eigen::Tensor<double, dim + 1>& get_global_u() { return global_u; }

        /**
         * @brief Get the global rho object
         * 
         * @return Eigen::Tensor<double, dim>& 
         */
        Eigen::Tensor<double, dim>& get_global_rho() { return global_rho; }

        /**
         * @brief Get the lattice dimensions object
         * 
         * @return std::array<Eigen::Index, dim> 
         */
        std::array<Eigen::Index, dim>& get_lattice_dimensions() { return lattice_dimensions; }

        /**
         * @brief Get the fluid nodes object
         * 
         * @return std::vector<Point<dim>>& 
         */
        std::vector<Point<dim>>& get_fluid_nodes(){ return fluid_nodes; }
        
        /**
         * @brief Get the boundary nodes object
         * 
         * @return std::vector<BoundaryPoint<dim>>& 
         */
        std::vector<BoundaryPoint<dim>>& get_boundary_nodes(){ return boundary_coord; }

        /**
         * @brief Get the inlet nodes object
         * 
         * @return std::vector<BoundaryPoint<dim>>& 
         */
        std::vector<BoundaryPoint<dim>>& get_inlet_nodes(){ return inlet_nodes_coord; }

        /**
         * @brief Get the outlet nodes object
         * 
         * @return std::vector<BoundaryPoint<dim>>& 
         */
        std::vector<BoundaryPoint<dim>>& get_outlet_nodes(){ return outlet_nodes_coord; }

        /**
         * @brief Get the obstacle nodes object
         * 
         * @return std::vector<BoundaryPoint<dim>>& 
         */
        std::vector<BoundaryPoint<dim>>& get_obstacle_nodes(){ return obstacle_nodes; }

        // ========================================================================================= 
        // ========================================================================================= 
        // ========================================================================================= 


    };
} // namespace llalbm::core

#endif // LLALBM_LATTICE_HPP