/**
 * @file MultiDimensionalLoop.hpp
 * @brief File containing a templated class that allows to loop over multiple dimensions, built around the Eigen::array object.
 */

#ifndef LLALBM_MULTIDIMENSIONALLOOP_HPP
#define LLALBM_MULTIDIMENSIONALLOOP_HPP

// =========== STL INCLUDES ===========
#include <iostream>
#include <type_traits>
// ====================================

// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Core>
#include <Eigen/Dense>
// ====================================

namespace llalbm::util
{

/**
 * @brief Base templated class for a multi-dimensional loop. It is specialized for the case in which the current dimension is the last one or not
 * using SFINAE.
 * 
 * @tparam ContainerObject Type of container that will store the grid positions, typically a vector.
 * @tparam DataObject Type of object that will represent each grid position. It needs a constructor that takes an Eigen::array.
 * @tparam dims Total number of dimensions.
 * @tparam current_dim Current dimension being processed in the loop.
 * @tparam Enabled SFINAE parameter to enable this specialization.
 */
template<typename ContainerObject, typename DataObject, std::size_t dims, std::size_t current_dim = 0, typename Enabled = void>
struct MultiDimensionalLoop;
 
/**
 * @brief Specialization of the MultiDimensionalLoop class for the case in which the current dimension is the last one. 
 * In this case, further recursive calls would reach a single element granularity, so a for loop is used to avoid 
 * the eventual overhead due to recursion
 * 
 * @tparam ContainerObject Type of container that will store the grid positions, typically a vector.
 * @tparam DataObject Type of object that will represent each grid position. It needs a constructor that takes an Eigen::array.
 * @tparam dims Total number of dimensions.
 * @tparam current_dim Current dimension being processed in the loop.
 * @tparam Enabled SFINAE parameter to enable this specialization.
 * 
 * @note Enabling of this specialization follows the SFINAE (Substitution Failure Is Not an Error) principle.
 */
template<typename ContainerObject, typename DataObject, std::size_t dims, std::size_t current_dim>
struct MultiDimensionalLoop<ContainerObject, DataObject, dims, current_dim, typename std::enable_if<current_dim == dims - 1>::type>
{
    static void assign_grid_positions(ContainerObject& container, Eigen::array<Eigen::Index, dims>& extensions, Eigen::array<Eigen::Index, dims>& grid_positions)
    {
        for (Eigen::Index i = 0; i < extensions[dims-1]; ++i)
        {
            grid_positions[dims-1] = i;
            container.push_back(DataObject(grid_positions));
        }
    }
};

/**
 * @brief A templated version of a multidimensional loop, used to access all the possible grid positions in a multi-dimensional grid.
 * Templating proves fundamental to allow scalability for a freely chosen number of dimensions.
 * 
 * @tparam ContainerObject Type of container that will store the grid positions, typically a vector.
 * @tparam DataObject Type of object that will represent each grid position. It needs a constructor that takes an Eigen::array.
 * @tparam dims Total number of dimensions.
 * @tparam current_dim Current dimension being processed in the loop.
 * @tparam Enabled SFINAE parameter to enable this specialization.
 * 
 */
template<typename ContainerObject, typename DataObject, std::size_t dims, std::size_t current_dim>
struct MultiDimensionalLoop <ContainerObject, DataObject, dims, current_dim, typename std::enable_if<current_dim < dims - 1>::type>
{
    static void assign_grid_positions(ContainerObject& container, Eigen::array<Eigen::Index, dims>& extensions, Eigen::array<Eigen::Index, dims>& grid_positions)
    {
        for (Eigen::Index i = 0; i < extensions[current_dim]; ++i)
        {
            grid_positions[current_dim] = i;
            MultiDimensionalLoop<ContainerObject, DataObject, dims, current_dim + 1>::assign_grid_positions(container, extensions, grid_positions);
        }
    }
};


}; // namespace llalbm::util

#endif // LLALBM_MULTIDIMENSIONALLOOP_HPP
