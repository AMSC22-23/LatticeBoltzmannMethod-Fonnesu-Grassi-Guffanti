#ifndef LLALBM_ALIASES_HPP
#define LLALBM_ALIASES_HPP

// =========== STL INCLUDES =============
#include <vector>
#include <array>
#include <cassert>
// ======================================

// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

using namespace Eigen; 


// =======================================================================================================
//              COORDINATES, NODES and BOUNDARIES
// =======================================================================================================

/// Generic coordinate of a point
using Coord = Eigen::Index;

/// @brief Different types of boundaries, identified when the lattice is created/read from file.
enum btypes {
    NONE = 0, 
    TOP_WALL,
    RIGHT_WALL,
    BOTTOM_WALL,
    LEFT_WALL,
    RIGHT_TOP_C,
    RIGHT_BOTTOM_C,
    LEFT_BOTTOM_C,
    LEFT_TOP_C
};

/// @brief Type of nodes
enum InputNodeType {
    FLUID = 0,
    SOLID = 1, 
    BOUNDARY = 2, 
    INLET = 3, 
    OUTLET = 4, 
    OBSTACLE = 5
};

/**
 * @brief A point in a d-dimensional space.
 * 
 * @tparam d number of dimensions of the domain.
 */
template<std::size_t d>
struct Point 
{
    /// @brief Coordinates of the point.
    Eigen::array<Eigen::Index, d> coords;

    /**
     * @brief Construct a new Point object
     * @note This is the default constructor.
     */
    Point() = default;

    /**
     * @brief Construct a new Point object by accepting an Eigen::array of coordinates.
     * 
     * @param coords_ array of coordinates.
     */
    Point(const Eigen::array<Eigen::Index, d>& coords_) : coords(coords_) {}

    /**
     * @brief Construct a new Point object by accepting an initializer list of coordinates
     * 
     * @param init_list 
     */
    Point(const std::initializer_list<Eigen::Index>& init_list)
    {
        assert(init_list.size() == d && "The number of coordinates must match the number of dimensions");
        std::copy(init_list.begin(), init_list.end(), coords.begin());
    }

    /**
     * @brief Operator overloading for comparison with another Point object.
     * 
     * @param other 
     * @return true if coordinates match
     * @return false if coordinates are not equal
     */
    bool operator==(const Point<d>& other) const
    {
        for(std::size_t i = 0; i < d; ++i)
        {
            if (coords[i] != other.coords[i])
                return false;
        }
        return true;
    }

    /**
     * @brief Operator overloading for comparison with another Point object.
     * 
     * @param other 
     * @return false if coordinates match
     * @return true if coordinates are not equal
     */
    bool operator!=(const Point<d>& other) const
    {
        return !(*this == other);
    }

    /**
    * @brief Operator overloading for comparison with another Point object.
    * The comparison, made by checking if the coordinates of the current object are less than the coordinates of the other object,
    * is used to sort the points in a container (more specifically the set into which points are placed when a lattice is being constructed in code).
    * 
    * @param other other point for comparison.
    */
    bool operator<(const Point<d>& other) const
    {
        for(std::size_t i = 0; i < d; ++i)
        {
            if (coords[i] < other.coords[i])
                return true;
            else if (coords[i] > other.coords[i])
                return false;
        }
        return false;
    }
};

/**
 * @brief Point at a boundary, with a type.
 * 
 * @tparam d number of dimensions of the domain.
 */
template<std::size_t d>
struct BoundaryPoint
{
    /// @brief Coordinates of the point.
    Eigen::array<Eigen::Index, d> coords;
    /// @brief  Type of boundary.
    btypes type;

    /**
     * @brief Construct a new Boundary Point object.
     * @note This is the default constructor.
     */
    BoundaryPoint() = default;

    /**
     * @brief Construct a new Boundary Point object accepting an array of coordinates.
     * 
     * @param coords_ array of coordinates.
     */
    BoundaryPoint(const Eigen::array<Eigen::Index, d>& coords_) : coords(coords_), type(btypes::NONE) {}

    /**
     * @brief Construct a new Boundary Point object accepting a Point object.
     * 
     * @param p point object whose coordinates are copied.
     */
    BoundaryPoint(const Point<d>& p) {
        for(std::size_t i = 0; i < d; ++i)
            coords[i] = p.coords[i];
        type = btypes::NONE;
    }
    /**
     * @brief Construct a new Boundary Point object accepting an initializer list of Eigen::Index values.
     * This method allows fast construction by passing initalizer lists, especially in tests.
     * 
     * @param init_list initializer list of Eigen::Index values.
     */
    BoundaryPoint(const std::initializer_list<Eigen::Index>& init_list)
    {
        assert(init_list.size() == d && "The number of coordinates must match the number of dimensions");
        std::copy(init_list.begin(), init_list.end(), coords.begin());
        type = btypes::NONE;
    }

    /**
     * @brief Operator overloading for comparison with another BoundaryPoint object.
     * 
     * @param other other boundary point object.
     * @return true if coordinates and type match.
     * @return false if coordinates or type are not equal.
     */
    bool operator==(const BoundaryPoint<d>& other) const
    {
        for(std::size_t i = 0; i < d; ++i)
        {
            if (coords[i] != other.coords[i])
                return false;
        }
        return type == other.type;
    }

    bool operator<(const BoundaryPoint<d>& other) const
    {
        for(std::size_t i = 0; i < d; ++i)
        {
            if (coords[i] < other.coords[i])
                return true;
            else if (coords[i] > other.coords[i])
                return false;
        }
        return false;
    }
};

/**
 * @brief Operator overloading for comparison between a Point and a BoundaryPoint object, independent of the type.
 * 
 * @tparam d number of dimensions of the domain.
 * @param p point object.
 * @param bp boundary point object.
 * @return true if coordinates match.
 * @return false if coordinates are not equal.
 */
template <std::size_t d>
bool operator==(const Point<d>& p, const BoundaryPoint<d>& bp)
{
    for(std::size_t i = 0; i < d; ++i)
    {
        if (p.coords[i] != bp.coords[i])
            return false;
    }
    return true;
}

// =======================================================================================================
// =======================================================================================================
// =======================================================================================================


// =======================================================================================================
//                  VELOCITY SETS (for 2D and 3D)
// =======================================================================================================

/// @brief D2Q9 velocity set
Matrix<double,9,3> D2Q9;


/// @brief D3Q19 velocity set
Matrix<double,19,4> D3Q19;

// =======================================================================================================
// =======================================================================================================
// =======================================================================================================


#endif