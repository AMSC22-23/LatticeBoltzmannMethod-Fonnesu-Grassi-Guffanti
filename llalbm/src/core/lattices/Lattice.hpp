#ifndef LLALBM_LATTICE_HPP
#define LLALBM_LATTICE_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
#include <cmath>
#include <map>
#include <memory>
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
#include "../../utils/analysis/FlowAnalyzer.hpp"
// =======================================


namespace llalbm::core
{
    using namespace llalbm::util::logger;
    using namespace llalbm::analysis;

    // Useful namespaces for PSBB
    using BlockCoord = std::pair<int, int>;
    using PointList = std::vector<ObstaclePoint<2>>;

    /**
     * The lattice class represents one of the most important components of a Lattice-Bolzmann simulation: the computational domain
     * onto which the simulation is performed. The class makes extensive use of the Eigen Linear Algebra library: the grid of points as
     * well as the particle distributions are stored as Eigen tensors, and are built based on the supplied spatial information when the object
     * is first instantiated. Templating allows for great versatility: with just a few limitations (mainly when using a GPU) the user is able to 
     * decide which policies to use or is free to implement other policies for the specific use-case.
     * 
     * @tparam ParallelizationPolicy Intermediary for calls to parallellized methods.
    */
    template<
        typename ParallelizationPolicy
    >
    class Lattice
    {
    private:
        static constexpr std::size_t dim = ParallelizationPolicy::dims;
        
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
        std::vector<ObstaclePoint<dim>> obstacle_nodes;

        /// @brief Number of velocities in the velocity set
        std::size_t q;

        // ========= POINTER TO THE FLOW ANALYZER =========
        std::shared_ptr<FlowAnalyzer<dim>> flow_analyzer = nullptr;

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
         * @brief Builds the tensors that describe the lattice.
         */
        void build_tensors()
        {
            // Print the dimensions
            logger.info("Dimensions:");
            for (std::size_t i = 0; i < dim; ++i) {
                logger.info("dim[" + std::to_string(i) + "] = " + std::to_string(lattice_dimensions[i]));
            }

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
         * @brief Construct a new completely empty lattice object
         * 
         */
        Lattice()
        :   q (0)
        ,   logger("LATTICE", std::cout)
        {}

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

            logger.info("Velocity Sets initialized");

            // Attach the nodes vectors to the initialization policy
            ParallelizationPolicy::attach_nodes(
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


            logger.info("Velocity Sets initialized");

            // Attach the nodes vectors to the initialization policy
            ParallelizationPolicy::attach_nodes(
                inlet_nodes_coord,
                outlet_nodes_coord
            );
            logger.info("Nodes attached to the initialization policy");

            // Then, lattice construction is complete.
            logger.info("Lattice is ready.");
        }

        /**
         * @brief Method called before starting the simulation when the partially saturated BB is used.
         * It computes weight parameters for each obstacle node.
         */
        void compute_obstacle_weight()
        {
            const std::size_t obstacles = obstacle_nodes.size();
            // Map that groups points in the same 5x5 square
            std::map<BlockCoord, PointList> blocks;

            for (std::size_t o_node = 0; o_node < obstacles; ++o_node)
            {
                const Eigen::Index i = obstacle_nodes[o_node].coords[0];
                const Eigen::Index j = obstacle_nodes[o_node].coords[1];

                // Define at which block the point is part of and push it in
                int block_x = i / 5;
                int block_y = j / 5;

                blocks[{block_x, block_y}].push_back(obstacle_nodes[o_node]);
            }

            // Let's calculate the weight b for each block
            for (auto &block : blocks)
            {
                const BlockCoord &block_coord = block.first;
                PointList &points_in_block = block.second;

                double num_elements = points_in_block.size();

                double b_value = num_elements / 25.0;

                for (auto &point : points_in_block)
                {
                    point.set_b(b_value);
                }
            }
        }   

        /**
         * @brief Executes the LBM simulation for the specified time.
         * 
         * @param time total time of the simulation
         * @param time_step time interval of the simulation
         * @param save_step number of iterations between each save
         * @param should_save whether Lattice should save the results
         */
        void perform_lbm(const double time, const double time_step = 1, const std::size_t save_step = 1, const bool should_save = true)
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
            
            if (flow_analyzer != nullptr)
            {
                flow_analyzer->initialize();
            }
            else
            {
                logger.info("Flow analyzer not attached");
            }
            
            // Compute the total number of iterations
            const std::size_t n_steps = time/time_step;
            std::size_t saved_file = 0;

            #ifdef LLALBM_VERBOSE
            logger.info("Initial equilibrium computation");
            #endif 

            // Inizitialization
            ParallelizationPolicy::calc_equilibrium(fluid_nodes, populations, global_u, global_rho);
            ParallelizationPolicy::calc_equilibrium(inlet_nodes_coord, populations, global_u, global_rho);
            ParallelizationPolicy::calc_equilibrium(outlet_nodes_coord, populations, global_u, global_rho);

            #ifdef LLALBM_VERBOSE
            logger.info("... Done!");
            #endif

            logger.info("Will do " + std::to_string(n_steps) + " steps");
            #ifdef TIMER
            auto start =std::chrono::high_resolution_clock::now();
            #endif
            for (std::size_t i = 0; i < n_steps; i++)
            {
                #ifndef LLALBM_VERBOSE
                std::cout << "\rStep " + std::to_string(i) + " / " + std::to_string(n_steps) << std::flush;
                #endif

                #ifdef LLALBM_VERBOSE
                logger.info("Step " + std::to_string(i) + " of " + std::to_string(n_steps));
                #endif 
                save = (i%save_step == 0 && should_save);


                // This will be controlled
                ParallelizationPolicy::update_nodes(i, global_u, global_rho);

                // 1) Update of the macroscopic quantities
                #ifdef LLALBM_VERBOSE
                logger.info("   Updating macroscopic quantities");
                #endif
                ParallelizationPolicy::update_macro(populations, fluid_nodes, global_rho, global_u);

                // 2) Compute equilibrium populations
                #ifdef LLALBM_VERBOSE
                logger.info("   Computing equilibrium populations");
                #endif
                ParallelizationPolicy::calc_equilibrium(fluid_nodes, equilibrium_populations,global_u,global_rho);
                ParallelizationPolicy::calc_equilibrium(inlet_nodes_coord, equilibrium_populations,global_u,global_rho);
                ParallelizationPolicy::calc_equilibrium(outlet_nodes_coord, equilibrium_populations,global_u,global_rho);


                
                // 3) Save ux and uy on files
                if(save)
                {
                    llalbm::util::writer::write_lattice_file<dim>(global_u, saved_file);
                    saved_file++;
                }

                // 4) Compute collisions
                #ifdef LLALBM_VERBOSE
                logger.info("   Colliding populations");
                #endif
                ParallelizationPolicy::collide(populations, equilibrium_populations, after_collision_populations, fluid_nodes, global_rho, global_u, time_step);
                ParallelizationPolicy::collide_open_inlet_boundary(populations, equilibrium_populations, after_collision_populations, inlet_nodes_coord, global_rho, global_u, time_step);
                ParallelizationPolicy::collide_open_outlet_boundary(populations, equilibrium_populations, after_collision_populations, outlet_nodes_coord, global_rho, global_u, time_step);

                if (flow_analyzer != nullptr && i > 0 && i % flow_analyzer->get_iterations_between_save() == 0)
                {
                    flow_analyzer->compute_flow_properties(after_collision_populations);
                }

                // 5) Propagate after collision populations, also to not fluid nodes
                #ifdef LLALBM_VERBOSE
                logger.info("   Streaming populations");
                #endif
                ParallelizationPolicy::stream(populations, after_collision_populations, fluid_nodes);
                ParallelizationPolicy::stream_open_inlet_boundary(populations, after_collision_populations, inlet_nodes_coord);
                ParallelizationPolicy::stream_open_outlet_boundary(populations, after_collision_populations, outlet_nodes_coord);


                // 6) Perform the collision at the boundaries
                #ifdef LLALBM_VERBOSE
                logger.info("   Computing Boundary action");
                #endif
                ParallelizationPolicy::update_domain_boundaries(populations, boundary_coord, global_rho, global_u);
                ParallelizationPolicy::update_obstacle_boundaries(populations, obstacle_nodes, global_rho, global_u);

                // 7) Perform the collision at the inlets and outlets
                #ifdef LLALBM_VERBOSE
                logger.info("   Computing Inlet and Outlet action");
                #endif
                ParallelizationPolicy::update_inlet_boundaries(populations, inlet_nodes_coord, global_rho, global_u);
                ParallelizationPolicy::update_outlet_boundaries(populations, outlet_nodes_coord, global_rho, global_u);

                #ifdef LLALBM_DEBUG
                std::cout << " -- Checking numerical stability" << std::endl;
                for (std::size_t i = 0; i < global_u.size(); i++)
                {
                    if (std::isnan(global_u.data()[i]))
                    {
                        logger.error("Numerical instability detected in the velocity tensor");
                        exit(1);
                    }
                }
                #endif

            }
            #ifdef TIMER
            auto end =std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end-start; 
            std::cout << std::endl << "Time: " << duration.count() << std::endl;
            #endif
            if (flow_analyzer != nullptr)
            {
                flow_analyzer->save_global_results();
            }


            logger.info("\n");
            logger.info("Done!");
        }

        /**
         * @brief Reinitializes a lattice object by reconstructing tensors and vectors of coordinates.
         * 
         * @param lattice_dimensions_ array indicating the extension of each spatial dimension in terms of grid points.
         * @param fluid_nodes_ list of coordinates of fluid nodes.
         * @param boundary_coord_ list of coordinates of boundary nodes.
         * @param inlet_nodes_coord_ list of coordinates of inlet nodes.
         * @param outlet_nodes_coord_ list of coordinates of outlet nodes.
         * @param obstacle_nodes_ list of coordinates of obstacle nodes.
         * 
         * @note This function is useful when the user wants to reinitialize the lattice with different data
         * or when the lattice has been created via the construction infrastructure.
         */
        void reinit(
            const std::size_t q_,
            const std::array<Eigen::Index, dim>& lattice_dimensions_,
            const std::vector<Point<dim>>& fluid_nodes_,
            const std::vector<BoundaryPoint<dim>>& boundary_coord_,
            const std::vector<BoundaryPoint<dim>>& inlet_nodes_coord_,
            const std::vector<BoundaryPoint<dim>>& outlet_nodes_coord_,
            const std::vector<ObstaclePoint<dim>>& obstacle_nodes_
        )
        {
            logger.info("Reinitializing lattice");
            q = q_;
            lattice_dimensions = lattice_dimensions_;
            build_tensors();

            boundary_coord = boundary_coord_;
            inlet_nodes_coord = inlet_nodes_coord_;
            outlet_nodes_coord = outlet_nodes_coord_;
            obstacle_nodes = obstacle_nodes_;
            fluid_nodes = fluid_nodes_;

            ParallelizationPolicy::attach_nodes(
                inlet_nodes_coord,
                outlet_nodes_coord
            );

            logger.info("Reinitialization completed");

        }

        // ========================================================================================= 
        //                                     GETTERS AND SETTERS     
        // ========================================================================================= 

        /**
         * @brief Get the populations object
         * 
         * @return Eigen::Tensor<double, dim + 1>& 
         */
        Eigen::Tensor<double, dim + 1>& get_populations() { return populations; }

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
         * @return std::vector<ObstaclePoint<dim>>& 
         */
        std::vector<ObstaclePoint<dim>>& get_obstacle_nodes(){ return obstacle_nodes; }

        /**
         * @brief Prints the lattice non-fluid structure to an output stream
         * 
         * @param out output stream onto which the structure is printed
         * @param print_fluid boolean indicating whether to print the fluid nodes coordinates
         */
        void print_lattice_structure(std::ostream& out, const bool print_fluid = false) const 
        {
            unsigned int n_nodes = 0;

            out << "Lattice Information" << std::endl;
            out << "-> List of obstacle nodes: " << std::endl;
            for (const auto& obstacle : obstacle_nodes)
            {
                for (unsigned int i = 0; i < dim; ++i)
                {
                    out << obstacle.coords[i] << " ";
                }
                out << "propagation direction " << obstacle.directions.to_string();
                out << std::endl;
                n_nodes++;
            }

            out << "-> List of inlet nodes: " << std::endl;
            for (const auto& inlet : inlet_nodes_coord)
            {
                for (unsigned int i = 0; i < dim; ++i)
                {
                    out << inlet.coords[i] << " ";
                }
                out << "type " << inlet.type;
                out << std::endl;
                n_nodes++;
            }

            out << "-> List of outlet nodes: " << std::endl;
            for (const auto& outlet : outlet_nodes_coord)
            {
                for (unsigned int i = 0; i < dim; ++i)
                {
                    out << outlet.coords[i] << " ";
                }
                out << "type " << outlet.type;
                out << std::endl;
                n_nodes++;
            }

            out << "-> List of boundary nodes: " << std::endl;
            for (const auto& boundary : boundary_coord)
            {
                for (unsigned int i = 0; i < dim; ++i)
                {
                    out << boundary.coords[i] << " ";
                }
                out << "type " << boundary.type;
                out << std::endl;
                n_nodes++;
            }

            out << "===================================" << std::endl;
            out << "Cumulative number of non fluid: " << n_nodes << std::endl;
            out << "===================================" << std::endl;

            if (print_fluid)
            {
                out << "List of fluid nodes: " << std::endl;
                for (const auto& fluid : fluid_nodes)
                {
                    for (unsigned int i = 0; i < dim; ++i)
                    {
                        out << fluid.coords[i] << " ";
                    }
                    out << std::endl;
                }
            
            }

        }

        /**
         * @brief Attaches a flow analyzer to the lattice object trhough a shared pointer
         * 
         * @param flow_analyzer_ std::shared_ptr<FlowAnalyzer<dim>>, instantiated by the user, to be attached to the lattice
         */
        void attach_flow_analyzer(std::shared_ptr<FlowAnalyzer<dim>> flow_analyzer_)
        {
            flow_analyzer = flow_analyzer_->shared_from_this();
        }
        // ========================================================================================= 
        // ========================================================================================= 
        // ========================================================================================= 


    };
} // namespace llalbm::core

#endif // LLALBM_LATTICE_HPP