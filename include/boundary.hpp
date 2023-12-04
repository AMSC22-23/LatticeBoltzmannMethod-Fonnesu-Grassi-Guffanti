#ifndef HH_BOUNDARY
#define HH_BOUNDARY

#include "latticeNode.hpp"
#include <vector>
#include <functional>


/**
 * Class that handles all possible methods to calculate the bouncing at the boundaries
 * @author Lorenzo Fonnesu
*/
class Boundary{
public:
    Boundary()=default;
    
    virtual ~Boundary() = default;
};

//We can use bounceBack when we have to deal with obstacle and 3D
class BounceBack : public Boundary
{
public:
    BounceBack() : Boundary(){};

    virtual ~BounceBack() = default;
};

/**
 * Implementation of the NEBB (Zou-He) method
 * @author Lorenzo Fonnesu
*/
class NEBB : public Boundary
{
public:
    NEBB() : Boundary(){};

    /**
     * Some useful constants
    */
    static constexpr double two_thirds = 2.0/3.0;
    static constexpr double one_sixth = 1.0/6.0;
    static constexpr double one_half = 0.5;

    /**
     * Function that handles the collision with a top wall
     * @param node whose data has to be updated
    */
    void calc_top_wall_bounce(LatticeNode<2> &node);

    /**
     * Function that handles the collision with a left wall
     * @param node whose data has to be updated
    */
    void calc_left_wall_bounce(LatticeNode<2> &node) const;

    /**
     * Function that handles the collision with a right wall
     * @param node whose data has to be updated
    */
    void calc_right_wall_bounce(LatticeNode<2> &node) const;

    /**
     * Function that handles the collision with a bottom wall
     * @param node whose data has to be updated
    */
    void calc_bottom_wall_bounce(LatticeNode<2> &node) const;

    /**
     * Function that handles the collision with a top left corner
     * @param node whose data has to be updated
    */
    void calc_top_left_corner_bounce(LatticeNode<2> &node) const;

    /**
     * Function that handles the collision with a top right corner
     * @param node whose data has to be updated
    */
    void calc_top_right_corner_bounce(LatticeNode<2> &node) const;

    /**
     * Function that handles the collision with a bottom left corner
     * @param node whose data has to be updated
    */
    void calc_bottom_left_corner_bounce(LatticeNode<2> &node) const;
     
    /**
     * Function that handles the collision with a right corner
     * @param node whose data has to be updated
    */
    void calc_bottom_right_corner_bounce(LatticeNode<2> &node) const;

    virtual ~NEBB() = default;
};

#endif