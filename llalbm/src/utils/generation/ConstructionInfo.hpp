/**
 * @file ConstructionInfo.hpp
 * @author Luca Guffanti
 * @brief Contains the definition and implementation of the ConstructionInfo class,
 * used to store the information needed to construct a Lattice object.
 */

#ifndef LLALBM_CONSTRUCTIONINFO_HPP
#define LLALBM_CONSTRUCTIONINFO_HPP

// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <algorithm>
#include <set>
// ======================================

// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "../../utils/loggers/Logger.hpp"
#include "../../utils/aliases.hpp" 
#include "../../utils/MultiDimensionalLoop.hpp"
// =======================================

namespace llalbm::util::generation
{
using namespace llalbm::util::logger;
using namespace llalbm::util;

enum NonFluidNodeType
{
    BOUNDARY = 0,
    INLET = 1,
    OUTLET = 2,
}; 


/**
 * @brief Class used to store information that will be used to construct a Lattice object.
 * This object is used to store intermediate information that will be used to directly construct a lattice,
 * that is, the dimensions of the computational domain, the list of inlet, outlet, boundary and obstacle nodes
 * expressed via their coordinates.
 * 
 * In order to acquire such information this object exposes a set of methods that allow the user to provide,
 * directly in code, intervals of coordinates that are associated to a specific note type.
 * When called, said functions make sure that no overlaps are present, that is that no point in the domain
 * is associated to more than one node type. Overlaps are solved by precedence, with the last caller being the
 * one that actually overwrites previous information, if present.
 * 
 * Last but not least, this object provides a method unload all the information, which is called by the Lattice constructor.
 * 
 * @tparam dim number of physical dimensions of the computational domain 
 */
template <std::size_t dim=2, std::size_t q=9>
class ConstructionInfo
{
private:

    /// Indicates if the dimensions of the computational domain have been provided.
    bool dimensions_provided = false;

    /// Indicates whether at least an inlet node has been provided.
    bool inlet_provided = false;

    /// Indicates whether at least an outlet node has been provided.
    bool outlet_provided = false;

    /// Indicates whether at least an obstacle node has been provided.
    bool obstacle_provided = false;

    /// Indicates whether lattice boundaries have been provided.
    bool boundaries_provided = false;

    /// List of dimensions of the computational domain. It's to be noted that coordinates the computational
    /// domain is modeled to be a hypercube, with edges extending from 0 to the dimensions.
    std::array<Eigen::Index, dim> domain_dimensions;

    /// List of boundary nodes
    std::set<BoundaryPoint<dim>> boundary_nodes;

    /// List of inlet nodes
    std::set<BoundaryPoint<dim>> inlet_nodes;

    /// List of outlet nodes
    std::set<BoundaryPoint<dim>> outlet_nodes;

    /// List of obstacle nodes
    std::set<ObstaclePoint<dim>> obstacle_nodes;

    /// Logger object
    Logger l;

public:

    ConstructionInfo(std::ostream& out_stream = std::cout)
    : l("ConstructionInfo", out_stream)
    {}

    /**
     * @brief Sets the extension of domain dimensions in the constructionInfo object
     * 
     * @param dims array of dimensions
     * @return number of elements that will constitute the lattice
     */
    Eigen::Index attach_domain_dimensions(const std::array<Eigen::Index, dim>& dims)
    {
        domain_dimensions = dims;
        
        // Make sure that no element of the dims array is equal to 0
        auto it = std::find(dims.begin(), dims.end(), 0);
        dimensions_provided = it == dims.end();
        const Eigen::Index elems = std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<Eigen::Index>())* dimensions_provided; 
        // And consequently return the total number of elements

        l.info("Provided domain with " + std::to_string(elems) + " elements");

        return elems;
    } 

    // ========================================================================================= 
    //                                INLES, OUTLETS AND BOUNDARIES
    //
    // The following methods are used to inlet, outlet or boundary nodes to the construction 
    // information used to build the lattice. Inlets, outlets and boundary are functionally
    // separated from obstacles, which are added with a different method as they must be 
    // treated differently.
    // ========================================================================================= 


    /**
     * @brief Adds an interval of nodes given the initial and final coordinates, whose type is chosen.
     * Nodes are added if and only if they are in a straight line and at the beginning or end (in any
     * dimension) of the physical domain.
     * The user is alerted in case of overlaps.
     * 
     * @return Eigen::Index number of added or modified nodes
     */
    Eigen::Index add_nodes_interval(const Point<dim>& begin, const Point<dim>& end, const NonFluidNodeType& type)
    {
        if (!dimensions_provided)
        {
            l.error("Dimensions of the computational domain have not been provided.");
            return 0;
        }

        if (!check_domain_bounds(begin) || !check_domain_bounds(end))
        {
            l.error("Provided points are not inside the domain");
            return 0;
        }

        std::size_t varying_dimension = 0;
        if (!are_in_line(begin, end, varying_dimension))
        {
            l.error("Provided points are not on the same line");
            return 0;
        }

        if (type == NonFluidNodeType::INLET)
        {
            inlet_provided = true;
        }
        else if (type == NonFluidNodeType::OUTLET)
        {
            outlet_provided = true;
        }
        else if (type == NonFluidNodeType::BOUNDARY)
        {
            boundaries_provided = true;
        }

        auto end_copy = end;
        end_copy.coords[varying_dimension]++;
        for (auto point = begin; point != end_copy; point.coords[varying_dimension]++)
        {
            if (check_possible_overlap<BoundaryPoint<dim>>(point, inlet_nodes))
            {
                #ifdef LLALBM_VERBOSE
                l.warn("Overlap detected in inlet nodes");
                #endif
                if (type == NonFluidNodeType::INLET)
                    continue;
                else
                    inlet_nodes.erase(point);
            }
            if (check_possible_overlap<BoundaryPoint<dim>>(point, boundary_nodes))
            {
                #ifdef LLALBM_VERBOSE
                l.warn("Overlap detected in inlet nodes");
                #endif
                if (type == NonFluidNodeType::BOUNDARY)
                    continue;
                else
                    boundary_nodes.erase(point);
            
            }
            if (check_possible_overlap<BoundaryPoint<dim>>(point, outlet_nodes))
            {
                #ifdef LLALBM_VERBOSE
                l.warn("Overlap detected in outlet nodes");
                #endif
                if (type == NonFluidNodeType::OUTLET)
                    continue;
                else
                    outlet_nodes.erase(point);

            }

            if (type == NonFluidNodeType::INLET)
                inlet_nodes.insert(BoundaryPoint<dim>(point));
            else if (type == NonFluidNodeType::OUTLET)
                outlet_nodes.insert(BoundaryPoint<dim>(point));
            else if (type == NonFluidNodeType::BOUNDARY)
                boundary_nodes.insert(BoundaryPoint<dim>(point));
        }
        

        const Eigen::Index added_nodes = std::abs(end.coords[varying_dimension] - begin.coords[varying_dimension]) + 1;
        l.info("Added " + std::to_string(added_nodes) + " nodes of type " + std::to_string(type));
        return added_nodes;

    }

    /**
     * @brief Constructs a perimeter of the specified type of node
     * 
     * @param type Type of node of the perimeter
     * @return Eigen::Index The number of added nodes
     */
    Eigen::Index add_perimeter_nodes(const NonFluidNodeType& type)
    {
        // First, check that the dimensions have been provided
        if (!dimensions_provided)
        {
            l.error("Dimensions of the computational domain have not been provided.");
            return 0;
        }

        // Then, add the perimeter nodes, which are points with one coordinate varying, and all the 
        // others being fixed to either 0 or the domain dimensions
        Eigen::Index added_nodes = 0;
        for (std::size_t i = 0; i < dim; i++)
        {
            Eigen::Index current_domain_dimension = domain_dimensions[i];
            
            // Starting with the coordinates set to 0
            for (Eigen::Index idx = 0; idx < current_domain_dimension; idx++)
            {
                std::array<Eigen::Index, dim> coords;
                coords.fill(0);
                coords[i] = idx;
                BoundaryPoint<dim> point(coords);
                insert_in_correct_set(point, type);
                added_nodes++;
            }

            std::array<Eigen::Index, dim> extensions;
            std::copy(domain_dimensions.begin(), domain_dimensions.end(), extensions.begin());
            std::for_each(extensions.begin(), extensions.end(), [](Eigen::Index& val){val = val-1;});
            // Then, with the coordinates set to the domain dimensions
            for (Eigen::Index idx = 0; idx < current_domain_dimension; idx++)
            {
                std::array<Eigen::Index, dim> coords;
                std::copy(extensions.begin(), extensions.end(), coords.begin());
                coords[i] = idx;
                BoundaryPoint<dim> point(coords);
                insert_in_correct_set(point, type);
                added_nodes++;
            }
        }
        return added_nodes;
    }


    // ========================================================================================= 
    //                                          OBSTACLES
    // Obstacles can be introduced in the lattice by using the following methods, which either 
    // introduce a shape (usually built by providing an origin point and an extension) or a
    // segment (built by providing two points) in the lattice.
    // 
    // ========================================================================================= 

    /**
     * @brief Adds an hyper rectangle to the computational domain, defined by an origin point and the extension of the solid along each dimension.
     * 
     * @param origin Origin point of the hyper rectangle. Considered to be the top left corner.
     * @param extensions Array of positive numbers indicating the extension of the hyper rectangle along each dimension.
     * @return Eigen::Index Number of added obstacle nodes.
     */
    Eigen::Index add_obstacle_hyper_rectangle(const Point<dim>& origin, const std::array<Eigen::Index,dim>& extensions)
    {
        // First verify that the domain has been constructed and that the point is inside the domain
        if (!dimensions_provided)
        {
            l.error("Dimensions of the computational domain have not been provided.");
            return 0;
        }

        if (!check_domain_bounds(origin))
        {
            l.error("Provided point is not inside the domain");
            return 0;
        }

        // Then check the other extreme of the hyper rectangle is within the domain
        Point<dim> end;
        for (std::size_t i = 0; i < dim; i++)
        {
            end.coords[i] = origin.coords[i] + extensions[i];
        }
        if (!check_obstacle_domain_bounds(end))
        {
            l.error("The hyper rectangle would get outside the domain");
            return 0;
        }

        // Then, add the obstacle nodes to the set of obstacle nodes
        Eigen::Index added_nodes = 0;
        std::vector<ObstaclePoint<dim>> obstacle_nodes_vec;


        Eigen::array<Eigen::Index, dim> grid_positions;
        Eigen::array<Eigen::Index, dim> extensions_eigen = extensions;
        Eigen::array<Eigen::Index, dim> origin_eigen = origin.coords;

        util::MultiDimensionalLoop<std::vector<ObstaclePoint<dim>>, ObstaclePoint<dim>, dim>::assign_grid_positions_with_offset(obstacle_nodes_vec, extensions_eigen, grid_positions, origin_eigen);
        std::move(obstacle_nodes_vec.begin(), obstacle_nodes_vec.end(), std::inserter(obstacle_nodes, obstacle_nodes.end()));
        // remove the all 0 point
        ObstaclePoint<dim> zero_point;

        zero_point.coords.fill(0);


        obstacle_nodes.erase(ObstaclePoint<dim>(zero_point));
        added_nodes = obstacle_nodes_vec.size();

        return added_nodes;
    }

    /**
     * @brief Adds a square obstacle to the computational domain, defined by an origin point and the extension of the solid along each dimension.
     * 
     * @param origin Origin point of the hyper square. Considered to be the top left corner.
     * @param extension Size of the square along a dimension.
     * @return Eigen::Index Number of added nodes.
     */
    Eigen::Index add_obstacle_hyper_square(const Point<dim>& origin, const Eigen::Index& extension)
    {
        
        // We can simply call the hyper rectangle generation with the extensions
        // being an array with all constants
        std::array<Eigen::Index, dim> ext;
        std::fill(ext.begin(), ext.end(), extension);

        return add_obstacle_hyper_rectangle(origin, ext);
    }

    /**
     * @brief Adds a sphere to the computational domain. The sphere is defined by its center and its radius.
     * 
     * @param origin Center of the sphere.
     * @param radius Radius of the sphere.
     * @return Eigen::Index Number of added nodes.
     */
    Eigen::Index add_obstacle_hyper_sphere(const Point<dim>& origin, const Eigen::Index& radius)
    {
        // First verify that the domain has been constructed and that the point is inside the domain
        if (!dimensions_provided)
        {
            l.error("Dimensions of the computational domain have not been provided.");
            return 0;
        }

        if (!check_domain_bounds(origin))
        {
            l.error("Provided point is not inside the domain");
            return 0;
        }

        // Then verify that the sphere would be built inside the domain.
        for (std::size_t i = 0; i < dim; i++) 
        {
            Point<dim> end_increased;
            Point<dim> end_decreased;
            std::copy(origin.coords.begin(), origin.coords.end(), end_increased.coords.begin());
            std::copy(origin.coords.begin(), origin.coords.end(), end_decreased.coords.begin());
            end_increased.coords[i] += radius;
            end_decreased.coords[i] -= radius;
            if (!check_obstacle_domain_bounds(end_increased) || !check_obstacle_domain_bounds(end_decreased))
            {
                l.error("The sphere would get outside the domain");
                return 0;
            }
        }

        // In order to add the obstacle nodes, the simplest thing is to build an hypersquare around the sphere,
        // and then remove the nodes that are outside the sphere by checking the euclidean distanc from the center.
        Eigen::Index added_nodes = 0;
        std::vector<ObstaclePoint<dim>> obstacle_nodes_vec;
        Eigen::array<Eigen::Index, dim> extensions;
        Eigen::array<Eigen::Index, dim> origin_eigen = origin.coords;
        Eigen::array<Eigen::Index, dim> grid_positions;

        std::fill(extensions.begin(), extensions.end(), radius);
        std::copy(origin.coords.begin(), origin.coords.end(), origin_eigen.begin());

        // But as the multidimensional loop works starting from the top left corner, we just need to translate the origin and adjust the extensions.
        std::transform(origin_eigen.begin(), origin_eigen.end(), extensions.begin(), origin_eigen.begin(), std::minus<Eigen::Index>());
        std::transform(extensions.begin(), extensions.end(), extensions.begin(), extensions.begin(), std::plus<Eigen::Index>());
        std::for_each(extensions.begin(), extensions.end(), [](Eigen::Index& val){val++;});


        util::MultiDimensionalLoop<std::vector<ObstaclePoint<dim>>, ObstaclePoint<dim>, dim>::assign_grid_positions_with_offset(obstacle_nodes_vec, extensions, grid_positions, origin_eigen);
        // For now simply add all the nodes to the set of obstacle nodes
        
        // Now we remove the nodes that are outside the sphere by iterating over the generated vector.
        const Eigen::Index radius_squared = radius*radius;
        std::vector<ObstaclePoint<dim>> to_remove;

        for (const auto& node : obstacle_nodes_vec)
        {  
            Eigen::array<Eigen::Index, dim> tmp;
            std::copy(node.coords.begin(), node.coords.end(), tmp.begin());
            std::transform(tmp.begin(), tmp.end(), origin.coords.begin(), tmp.begin(), std::minus<Eigen::Index>());
            const std::size_t distance_squared = std::inner_product(tmp.begin(), tmp.end(), tmp.begin(), 0);
            if (distance_squared > radius_squared)
            {
                to_remove.push_back(node);
            }
        }

        // And finally remove the elements that are outside the sphere
        for (const auto& node : to_remove)
        {
            obstacle_nodes_vec.erase(std::remove(obstacle_nodes_vec.begin(), obstacle_nodes_vec.end(), node), obstacle_nodes_vec.end());
        }


        added_nodes = obstacle_nodes_vec.size();
        std::move(obstacle_nodes_vec.begin(), obstacle_nodes_vec.end(), std::inserter(obstacle_nodes, obstacle_nodes.end()));
        return added_nodes;
    }

    
    // ========================================================================================= 
    //                                      GETTERS AND SETTERS
    // ========================================================================================= 

    /**
     * @brief Returns the dimensions of the computational domain.
     * 
     * @return std::array<Eigen::Index, dim> 
     */
    std::array<Eigen::Index, dim> get_domain_dimensions() const
    {
        return domain_dimensions;
    }

    /**
     * @brief Returns the list of boundary nodes.
     * 
     * @return std::set<BoundaryPoint<dim>> 
     */
    std::set<BoundaryPoint<dim>> get_boundary_nodes() const
    {
        return boundary_nodes;
    }

    /**
     * @brief Returns the list of inlet nodes.
     * 
     * @return std::set<BoundaryPoint<dim>> 
     */
    std::set<BoundaryPoint<dim>> get_inlet_nodes() const
    {
        return inlet_nodes;
    }

    /**
     * @brief Returns the list of outlet nodes.
     * 
     * @return std::set<BoundaryPoint<dim>> 
     */
    std::set<BoundaryPoint<dim>> get_outlet_nodes() const
    {
        return outlet_nodes;
    }

    /**
     * @brief Returns the list of obstacle nodes.
     * 
     * @return std::set<BoundaryPoint<dim>> 
     */
    std::set<ObstaclePoint<dim>> get_obstacle_nodes() const
    {
        return obstacle_nodes;
    }

    bool are_dimensions_provided() const
    {
        return dimensions_provided;
    }

    bool are_inlets_provided() const
    {
        return inlet_provided;
    }

    bool are_outlets_provided() const
    {
        return outlet_provided;
    }

    bool are_obstacles_provided() const
    {
        return obstacle_provided;
    }

    bool are_boundaries_provided() const
    {
        return boundaries_provided;
    }

    // ========================================================================================= 
    // ========================================================================================= 
    // =========================================================================================

private:

    // TODO: WHEN ADDING CLOSED SHAPE OBSTACLES, THE INSIDE MUST BE SET TO BE AN OBSTACLE


    /**
     * @brief Checks whether two points are in a straight line
     *  
     * @returns True of false based on the check
     */
    bool are_in_line(const Point<dim> first, const Point<dim> second, std::size_t& varying_dimension)
    {
        bool possible_varying = false;
        for (unsigned int i = 0; i < dim; i++)
        {
            if (first.coords[i] != second.coords[i])
            {
                if (possible_varying)
                {
                    l.error("Provided points are not on the same line");
                    return false;
                }
                possible_varying = true;
                varying_dimension = i;
            }
        }
        return true;        
    }

    /**
     * @brief Checks whether the provided coordinate are within the computational domain.
     * 
     * @returns True or false whether based on the check.
     */
    bool check_domain_bounds(const Point<dim>& point) const
    {
        for (unsigned int i = 0; i < dim; i++)
        {
            if(point.coords[i] < 0 || point.coords[i] >= domain_dimensions[i])
            {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Checks whether the provided coordinate are within the computational domain.
     * 
     * @returns True or false whether based on the check.
     */
    bool check_obstacle_domain_bounds(const Point<dim>& point) const
    {
        for (unsigned int i = 0; i < dim; i++)
        {
            if(point.coords[i] <= 0 || point.coords[i] >= domain_dimensions[i] - 1)
            {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Checks whether the provided coordinates overlap with others that are registered as non-generic fluid nodes. Even though this method
     * is not necessarily vital for the correct construction because overlaps are solved by picking the last caller, it's still useful 
     * to log the presence of overlaps.
     * 
     * @param point Coordinates of the point for which the check is performed.
     * @param to_check Set of coordinates for which the check must be done
     * 
     * @returns Whether an ovelap is present.
     */
    template<typename ObjectType>
    bool check_possible_overlap(const ObjectType& point, const std::set<ObjectType>& to_check) const
    {
        const auto it = to_check.find(point);
        if (it != to_check.end())
        {
            return true;
        }
        return false;
    }

    /**
     * @brief Inserts a @link BoundaryPoint @endlink in the correct set, based on the type of node.
     * 
     * @param point Point to be inserted in the set
     * @param type Type of the point
     */
    void insert_in_correct_set(const BoundaryPoint<dim>& point, const NonFluidNodeType& type)
    {
        if (type == NonFluidNodeType::INLET)
        {
            inlet_nodes.insert(point);
        }
        else if (type == NonFluidNodeType::OUTLET)
        {
            outlet_nodes.insert(point);
        }
        else if (type == NonFluidNodeType::BOUNDARY)
        {
            boundary_nodes.insert(point);
        }
    }
};
}; // namespace llalbm::util::generation
#endif // LLALBM_CONSTRUCTIONINFO_HPP