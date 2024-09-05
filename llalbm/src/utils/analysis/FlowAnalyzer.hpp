/**
 * @file FlowAnalyzer.hpp
 * @author Luca GUffanti
 * @brief Contains the FlowAnalyzer class, used to compute the total lift of the obstacle nodes.
 * 
 * 
 * 
 * 
 */

#ifndef LLALBM_FLOW_ANALYZER_HPP
#define LLALBM_FLOW_ANALYZER_HPP

// ======= LLALBM INCLUDES =======
#include "../aliases.hpp"
#include "../MultiDimensionalLoop.hpp"
#include "../loggers/Logger.hpp"
// ===============================

// ======= STL INCLUDES ==========
#include <vector>
#include <array>
#include <unordered_set>
#include <memory>
#include <filesystem>
// ===============================

// ======= EIGEN INCLUDES ========
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
// ===============================


namespace llalbm::analysis
{

using namespace llalbm::util::logger;

/**
 * @brief Class used to compute the lift and drag acting on the obstacle nodes, using the momentum exchange approach.
 * 
 * The class is templated on the dimension of the problem, with specialization for 2 dimensions.
 * Specialization is required because the way lift and drag are computed depends on the choice of velocity set for the problem,
 * so the non-specialized version of the class is not implemented.
 * 
 * @tparam dim dimension of the problem
 */
template <std::size_t dim>
class FlowAnalyzer : public std::enable_shared_from_this<FlowAnalyzer<dim>>
{
private:
    Logger l;
    double total_lift;
    double total_drag;
    std::size_t iterations_between_save;
    int save_counter;

public:
    void add_point(const ObstaclePoint<dim>& point)
    {
        throw std::runtime_error("FlowAnalyzer is not implemented for this dimension");
    }

    void add_point_vec(std::vector<ObstaclePoint<dim>>& point_vec)
    {
        throw std::runtime_error("FlowAnalyzer is not implemented for this dimension");
    }

    void clear_points()
    {
        throw std::runtime_error("FlowAnalyzer is not implemented for this dimension");
    }

    void recognize_nearby_obstacle(const ObstaclePoint<dim>& origin, std::vector<ObstaclePoint<dim>>& obstacle_points)
    {
        throw std::runtime_error("FlowAnalyzer is not implemented for this dimension");
    }

    std::pair<double,double> compute_flow_properties(const Eigen::Tensor<double, dim+1>& densities, const bool& should_save=false)
    {
        throw std::runtime_error("FlowAnalyzer is not implemented for this dimension");
    }

    void set_iterations_between_save(const std::size_t& iterations_between_save_)
    {
        iterations_between_save = iterations_between_save_;
    }

    std::size_t get_iterations_between_save() const
    {
        return iterations_between_save;
    }

    void initialize() 
    {
        throw std::runtime_error("FlowAnalyzer is not implemented for this dimension");
    }

    std::shared_ptr<FlowAnalyzer<dim>> create()
    {
        return nullptr;
    }
};

/**
 * @brief Template specialization for the 2D case. 
 * The FlowAnalyzer manages the computation, saving and output of the lift and drag forces exerted on obstacle nodes.
 * Such obstacles are passed by the user directly to the FlowAnalyzer object, which can be used in two ways:
 * 
 * 1. To perform analysis at the end of simulations: in this case it is enough to call the compute_flow_properties method
 * as the simulation ends passing the population tensor. The method will return the total lift and drag forces, eventually saving
 * the results to ouptut files.
 * 
 * 2. To perform analysis at regular intervals during the simulation: in this case the FlowAnalyzer must be attached to the lattice object 
 * by wrapping it inside a shared pointer. The lattice object will then take care of computation and saving of the results.
 * 
 * It is to be noted that there are also many ways of passing obstacle nodes to the LiftAnalyzer:
 * 
 * 1. By adding them one by one using the add_point method
 * 
 * 2. By adding a list of them using the add_point_vec method
 * 
 * 3. By allowing the FlowAnalyzer to automatically compute the obstacle nodes linked to a given origin point using the recognize_nearby_obstacle method.
 * In this last case, it is then fundamental to call the match_wet_nodes method to filter out the obstacle nodes that are not linked to fluid nodes, if the
 * list of nodes passed to recognize_nearby_obstacle is not already filtered.
 * 
 */
template<>
class FlowAnalyzer<2> : public std::enable_shared_from_this<FlowAnalyzer<2>>
{
private:
    /**
     * @brief List of all the obstacle points that are considered for the analysis.
     * Such obstacle points are linked to liquid nodes, to which populations are propagated using 
     * a bounce back approach.
     * 
     */
    std::vector<ObstaclePoint<2>> considered_points;

    /**
     * @brief Logger object used to print messages
     * 
     */
    Logger l;

    /**
     * @brief Total lift, computed as the sum of the vertical forces acting on the obstacle nodes.
     * 
     */
    double total_lift;

    /**
     * @brief Total drag, computed as the sum of the horizontal forces acting on the obstacle nodes.
     * 
     */
    double total_drag;

    /**
     * @brief How many iterations should pass between each save of the contribution of each node to lift and drag.
     * 
     */
    std::size_t iterations_between_save;

    /**
     * @brief How many iterations have been saved so far.
     * 
     */
    int save_counter;


    /**
     * @brief Whether or not data from the analysis should be saved, default is true.
     * 
     */
    bool should_save = true;

    /**
     * @brief Output file stream object for total lift and drag contributions.
     * 
     */
    std::ofstream total_lift_drag_out;

    /**
     * @brief Vector containing pairs of total lift and drag for each iteration.
     */
    std::vector<std::pair<double,double>> lift_drag_per_iteration;

public: 


    FlowAnalyzer(const std::size_t& iterations_between_save_=10) :
    l ("FlowAnalyzer", std::cout),
    total_lift(0.),
    total_drag(0.),
    iterations_between_save(iterations_between_save_)
    {
        initialize();
    }

    FlowAnalyzer(std::vector<ObstaclePoint<2>>& considered_points_, const std::size_t& iterations_between_save_=10, const bool& should_save_=true) :
    l ("FlowAnalyzer", std::cout),
    considered_points(considered_points_),
    total_lift(0.),
    total_drag(0.),
    should_save(should_save_),
    iterations_between_save(iterations_between_save_)
    {
        initialize();
    }

    /**
     * @brief Factory method to create a shared pointer to a FlowAnalyzer object.
     * 
     * @param considered_points_ considered points vector (must not be temporary)
     * @param iterations_between_save_ iterations between save value.
     * @return std::shared_ptr<FlowAnalyzer<2>> 
     */
    static std::shared_ptr<FlowAnalyzer<2>> create(std::vector<ObstaclePoint<2>>& considered_points_,
        const std::size_t& iterations_between_save_,
        const bool& should_save_) 
    {
        return std::make_shared<FlowAnalyzer<2>>(considered_points_, iterations_between_save_, should_save_);
    }

    /**
     * @brief Initializes the flow analyzer by constructing the empty output directory.
     * The initialization function is performed only once, at the beginning of the simulation, and 
     * has been separated from the constructor to allow for more flexibility (i.e. direct calling from the lattice object).
     */
    void initialize()
    {
        l.info("Initializing FlowAnalyzer");
        save_counter = 0;
        if (!std::filesystem::create_directory("./flow_analysis"))
        {
            l.info("Directory already exists");
            // If the directory already exists, clear it
            std::vector<std::filesystem::path> to_remove;
            for (const auto& entry : std::filesystem::directory_iterator("./flow_analysis"))
            {
                if (entry.path().extension() == ".csv")
                {
                    to_remove.push_back(entry.path());
                }
            } 

            for (const auto& path : to_remove)
            {
                std::filesystem::remove(path);
            }
        }
    }

    /**
     * @brief Saves the global lift and drag results vector in a csv file. This approach has been selected to avoid keeping
     * the output file open during the entirety of the simulation, as this affects performance.
     * 
     */
    void save_global_results() const 
    {
        std::ofstream total_lift_drag_out;
        total_lift_drag_out.open("./flow_analysis/total_lift_drag.csv");
        total_lift_drag_out << "iteration,total_lift,total_drag" << std::endl;
        for (std::size_t i = 0; i < lift_drag_per_iteration.size(); ++i)
        {
            total_lift_drag_out << i << "," << lift_drag_per_iteration[i].first << "," << lift_drag_per_iteration[i].second << std::endl;
        }
    }

    /**
     * @brief Add a point to the list of obstacle points that are considered for the analysis.
     * 
     * @param point obstacle point to be added
     */
    void add_point(const ObstaclePoint<2>& point)
    {
        considered_points.push_back(point);
    }

    /**
     * @brief Add a vector of points to the list of obstacle points that are considered for the analysis.
     * 
     * @param point_vec vector of obstacle points to be added
     */
    void add_point_vec(std::vector<ObstaclePoint<2>>& point_vec)
    {
        considered_points = point_vec;
    }

    /**
     * @brief Clears the vector of obstacle points that are considered for the analysis.
     * 
     */
    void clear_points()
    {
        considered_points.clear();
    }

    /**
     * @brief Set the iterations between save to file.
     * 
     * @param iterations_between_save_ How many lattice boltzmann iterations should pass between each save of the contributions.
     */
    void set_iterations_between_save(const std::size_t& iterations_between_save_)
    {
        iterations_between_save = iterations_between_save_;
    }

    void set_should_save(const bool& should_save_)
    {
        should_save = should_save_;
    }

    std::size_t get_iterations_between_save() const
    {
        return iterations_between_save;
    }

    /**
     * @brief Given an origin point, this function identifies all obstacle points that are linked to the origin point.
     * Two nodes are linked if there exists a path of nodes between the two.
     * 
     * @param origin origin point for the analysis
     * @param obstacle_points list of obstacle points that are checked to the origin point
     */
    void recognize_nearby_obstacle(const ObstaclePoint<2>& origin, std::vector<ObstaclePoint<2>>& obstacle_points)
    {
        std::vector<ObstaclePoint<2>> obstacle_points_copy(obstacle_points.size());
        std::copy(obstacle_points.begin(), obstacle_points.end(), obstacle_points_copy.begin());
        l.info("Computing obstacles");
        compute_nearby_obstacles(origin, obstacle_points_copy);
    }


    /**
     * @brief Computes the flow properties, that are lift and drag. This function should be called after the collision but before
     * the streaming, as computation simulates the momentum exchange in the case of a bounce back policy
     * using the momentum exchange approach. Eventually, single contributions to lift and drag can be saved in a csv file.
     *  
     * Nonetheless, input densities are not modified by this function, so it can be called at any time.
     * 
     * Lift is the total vertical force acting on the obstacle, while drag is the total horizontal force acting on the obstacle.
     * 
     * @param densities populations tensor
     * 
     * @returns a pair containing the total lift and the total drag
     */
    std::pair<double,double> compute_flow_properties(const Eigen::Tensor<double, 3>& densities)
    {
#ifdef LLALBM_VERBOSE
        l.info("Computing flow properties");
#endif 
        std::ofstream out_stream;

        if (should_save)
        {
#ifdef LLALBM_VERBOSE
            l.info("Saving contributions to lift and drag");
#endif
            std::string filename = "./flow_analysis/contribution_"+std::to_string(save_counter)+".csv";
            out_stream.open(filename);
            out_stream << "i,j,lift,drag" << std::endl;
        }

        total_lift = 0.;
        total_drag = 0.;

        for(const auto& node : considered_points)
        {
            double local_lift = 0.;
            double local_drag = 0.;
            const Eigen::Index i = node.coords[0];
            const Eigen::Index j = node.coords[1];

            // Based on the propagation direction, we can understand where fluid nodes are present.
            //               6 2 5
            //                \|/
            //              3--0--1
            //                /|\
            //               7 4 8
            // Where this is centered on the obstacle node. This means that
            // -> If 1, then propagation on E, and momentum is applied from W 
            // -> If 2, then propagation above, so the fluid node is above (North) and momentum is applied from the top through 4
            // -> If 3, then propagation on W, and momentum is applied from E 
            // -> If 4, then propagation below, so the fluid node is below (South) and momentum is applied from below through 2
            // -> If 5, then propagation on NE, and momentum is applied from SW through 7
            // -> If 6, then propagation on NW, and momentum is applied from SE through 8
            // -> If 7, then propagation on SW, and momentum is applied from NE through 5
            // -> If 8, then propagation on SE, and momentum is applied from NW through 6
            // All in all:
            //  Contribution to LIFT : populations 2, 4, 5, 6, 7, 8
            //  Contribution to DRAG : populations 1, 3, 5, 6, 7, 8

            if (node.directions.test(1))
            {
                local_drag -= densities(i, j+1, 3);
            }
            if (node.directions.test(2)) 
            {
                local_lift -= densities(i-1, j, 4);
            }
            if (node.directions.test(3))
            {
                local_drag += densities(i, j-1, 1);
            }
            if (node.directions.test(4)) 
            {
                local_lift += densities(i+1, j, 2);
            }
            if (node.directions.test(5)) 
            {
                local_drag -= (densities(i-1, j+1, 7)) * std::sqrt(2)/2.0;
                local_lift -= (densities(i-1, j+1, 7)) * std::sqrt(2)/2.0;

            }
            if (node.directions.test(6)) 
            {
                local_drag += (densities(i-1, j-1, 8)) * std::sqrt(2)/2.0;
                local_lift -= (densities(i-1, j-1, 8)) * std::sqrt(2)/2.0;   
            }
            if (node.directions.test(7)) 
            {
                local_drag += (densities(i+1, j-1, 5)) * std::sqrt(2)/2.0;
                local_lift += (densities(i+1, j-1, 5)) * std::sqrt(2)/2.0;
            }
            if (node.directions.test(8)) 
            {
                local_drag -= (densities(i+1, j+1, 6)) * std::sqrt(2)/2.0;
                local_lift += (densities(i+1, j+1, 6)) * std::sqrt(2)/2.0;
            }

            if (should_save)
            {
                out_stream << i << "," << j << "," << local_lift << "," << local_drag << std::endl;
            }

            total_drag += local_drag;
            total_lift += local_lift;
        }

        if (should_save)
        {
            total_lift_drag_out << save_counter << "," << total_lift << "," << total_drag << std::endl;
#ifdef LLALBM_VERBOSE
            l.info("Contributions to lift and drag saved in contribution.csv");
#endif
            save_counter++;
        }
        lift_drag_per_iteration.push_back(std::make_pair(total_lift, total_drag));
        return std::make_pair(total_lift, total_drag);
    }

    /**
     * @brief This function modifies the considered points, by only selecting those that are wet, that is, 
     * that would be propagating populations to the fluid nodes. In this way, the lift and drag are computed only
     * on useful nodes.
     * 
     * @param wet_obstacle_nodes set of obstacle nodes that propagate to fluid nodes.
     */
    void match_wet_nodes(const std::vector<ObstaclePoint<2>>& wet_obstacle_nodes)
    {
        std::vector<ObstaclePoint<2>> new_nodes;
        // Iterate over all considered points and check if they are in the wet_obstacle_nodes.
        // If they are, then add them to the new_nodes list
        for (auto& it: considered_points)
        {
            if (std::find(wet_obstacle_nodes.begin(), wet_obstacle_nodes.end(), it) != wet_obstacle_nodes.end())
            {
                new_nodes.push_back(it);
            }
        }

        // Update the considered points
        considered_points = new_nodes;  
    }

    /**
     * @brief Getter for the considered points
     * 
     * @return std::vector<ObstaclePoint<2>> considered points
     */
    std::vector<ObstaclePoint<2>> get_points() const
    {
        return considered_points;
    }

    /**
     * @brief Get the total lift 
     * 
     * @return double 
     */
    double get_total_lift() const
    {
        return total_lift;
    }

    /**
     * @brief Get the total drag
     * 
     * @return double 
     */
    double get_total_drag() const
    {
        return total_drag;
    }

    /**
     * @brief Print the considered points for analysis to the default output stream or to a csv file.
     * 
     * @param print_to_file Flag indicating whether to print to a file or to the default output stream
     */
    void print_considered_points(const bool& print_to_file=true) const
    {
        if (!print_to_file)
        {
            for (const auto& node : considered_points)
            {
                std::cout << node.coords[0] << "," << node.coords[1] << std::endl;
            }
        }
        else
        {
            std::ofstream out("obstacle_layout.csv");
            out << "x,y" << std::endl;
            for (const auto& node : considered_points)
            {
                out << node.coords[0] << "," << node.coords[1] << std::endl;
            }
            l.info("Obstacle layout saved in obstacle_layout.csv");
        }
    }

    /**
     * @brief Get all the global lift and drag results for each iteration.
     * 
     * @return std::vector<std::pair<double,double>> 
     */
    std::vector<std::pair<double,double>> get_all_results() const
    {
        return lift_drag_per_iteration;
    }

private: 
    /**
     * @brief Recursive method used to compute all neighboring obstacles, starting from a given origin point.
     * @note This method is kickstarted by the recognize_nearby_obstacles method.
     * 
     * @param origin origin point
     * @param obstacle_points_copy copy of the obstacle points list 
     */
    void compute_nearby_obstacles(const ObstaclePoint<2>& origin, std::vector<ObstaclePoint<2>>& obstacle_points_copy)
    {
        // Find the origin in the list and add it to the considered points if not already present
        auto it = std::find(obstacle_points_copy.begin(), obstacle_points_copy.end(), origin);
        if (it != obstacle_points_copy.end())
        {
            considered_points.push_back(origin);
            obstacle_points_copy.erase(it);
        }

        // Iterate over all neighboring cells (including diagonals)
        for (int i = -1; i <= 1; ++i)
        {
            for (int j = -1; j <= 1; ++j)
            {
                // Skip the origin itself
                if (i == 0 && j == 0) continue;

                ObstaclePoint<2> neighbor = origin;
                neighbor.coords[0] += i;
                neighbor.coords[1] += j;

                // Check if this neighbor is in the obstacle_points_copy list
                auto neighbor_it = std::find(obstacle_points_copy.begin(), obstacle_points_copy.end(), neighbor);
                if (neighbor_it != obstacle_points_copy.end())
                {
                    considered_points.push_back(neighbor);
                    obstacle_points_copy.erase(neighbor_it);
                    compute_nearby_obstacles(neighbor, obstacle_points_copy);
                }
            }
        }
    }

};

}

#endif //LLALBM_FLOWANALYZER_HPP