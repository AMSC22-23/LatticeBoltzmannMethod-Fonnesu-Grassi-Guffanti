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
    OBSTACLE = 3
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
template <unsigned int dim>
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
    Eigen::array<Eigen::Index, dim> domain_dimensions;

    /// List of boundary nodes
    std::set<Point<dim>> boundary_nodes;

    /// List of inlet nodes
    std::set<Point<dim>> inlet_nodes;

    /// List of outlet nodes
    std::set<Point<dim>> outlet_nodes;

    /// List of obstacle nodes
    std::set<Point<dim>> obstacle_nodes;

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
    Eigen::Index attach_domain_dimensions(const Eigen::array<Eigen::Index, dim>& dims)
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
    //                                     INLES AND OUTLETS
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
        else if (type == NonFluidNodeType::OBSTACLE)
        {
            obstacle_provided = true;
        }
        else if (type == NonFluidNodeType::BOUNDARY)
        {
            boundaries_provided = true;
        }

        // Check for overlaps: the point must not be present in any other list
        for (auto point = begin; point != end; point.coords[varying_dimension]++)
        {
            if (check_possible_overlap(point, inlet_nodes))
            {
                l.warn("Overlap detected in inlet nodes");
                if (type == NonFluidNodeType::INLET)
                    continue;
                else
                    inlet_nodes.erase(point);
            }
            if (check_possible_overlap(point, boundary_nodes))
            {
                l.warn("Overlap detected in boundary nodes");
                if (type == NonFluidNodeType::BOUNDARY)
                    continue;
                else
                    boundary_nodes.erase(point);
            
            }
            if (check_possible_overlap(point, outlet_nodes))
            {
                l.warn("Overlap detected in outlet nodes");
                if (type == NonFluidNodeType::OUTLET)
                    continue;
                else
                    outlet_nodes.erase(point);

            }
            if (check_possible_overlap(point, obstacle_nodes))
            {
                l.warn("Overlap detected in obstacle nodes");
                if (type == NonFluidNodeType::OBSTACLE)
                    continue;
                else
                    obstacle_nodes.erase(point);

            }
            

            if (type == NonFluidNodeType::INLET)
                inlet_nodes.insert(point);
            else if (type == NonFluidNodeType::OUTLET)
                outlet_nodes.insert(point);
            else if (type == NonFluidNodeType::OBSTACLE)
                obstacle_nodes.insert(point);
            else if (type == NonFluidNodeType::BOUNDARY)
                boundary_nodes.insert(point);
        }
        

        const Eigen::Index added_nodes = std::abs(end.coords[varying_dimension] - begin.coords[varying_dimension]) + 1;
        l.info("Added " + std::to_string(added_nodes) + " nodes of type " + std::to_string(type));
        return added_nodes;

    }

    
    // ========================================================================================= 
    //                                      GETTERS AND SETTERS
    // ========================================================================================= 

    /**
     * @brief Returns the dimensions of the computational domain.
     * 
     * @return Eigen::array<Eigen::Index, dim> 
     */
    Eigen::array<Eigen::Index, dim> get_domain_dimensions() const
    {
        return domain_dimensions;
    }

    /**
     * @brief Returns the list of boundary nodes.
     * 
     * @return std::set<Point<dim>> 
     */
    std::set<Point<dim>> get_boundary_nodes() const
    {
        return boundary_nodes;
    }

    /**
     * @brief Returns the list of inlet nodes.
     * 
     * @return std::set<Point<dim>> 
     */
    std::set<Point<dim>> get_inlet_nodes() const
    {
        return inlet_nodes;
    }

    /**
     * @brief Returns the list of outlet nodes.
     * 
     * @return std::set<Point<dim>> 
     */
    std::set<Point<dim>> get_outlet_nodes() const
    {
        return outlet_nodes;
    }

    /**
     * @brief Returns the list of obstacle nodes.
     * 
     * @return std::set<Point<dim>> 
     */
    std::set<Point<dim>> get_obstacle_nodes() const
    {
        return obstacle_nodes;
    }

    // ========================================================================================= 
    // ========================================================================================= 
    // =========================================================================================

private:

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
     * @brief Checks whether the provided coordinates overlap with others that are registered as non-generic fluid nodes. Even though this method
     * is not necessarily vital for the correct construction because overlaps are solved by picking the last caller, it's still useful 
     * to log the presence of overlaps.
     * 
     * @param point Coordinates of the point for which the check is performed.
     * @param to_check Set of coordinates for which the check must be done
     * 
     * @returns Whether an ovelap is present.
     */
    bool check_possible_overlap(const Point<dim>& point, const std::set<Point<dim>>& to_check) const
    {
        const auto it = to_check.find(point);
        if (it != to_check.end())
        {
            return true;
        }
        return false;
    }
};

};

#endif // LLALBM_CONSTRUCTIONINFO_HPP