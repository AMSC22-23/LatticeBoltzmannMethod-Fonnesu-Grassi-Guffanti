#ifndef HH_BOUNDARY
#define HH_BOUNDARY

#include <vector>
#include <functional>


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
     * @param f vector containing distribution functions of a node
     * @param uy0 vertical velocity of the node
     * @return the pre-streaming vector of distribution functions of the node 
    */
    std::vector<double> calc_top_wall_bounce(std::vector<double> &f, double uy0) const;

    /**
     * @param f vector containing distribution functions of a node
     * @param uy0 vertical velocity of the node
     * @return the pre-streaming vector of distribution functions of the node 
    */
    std::vector<double> calc_left_wall_bounce(std::vector<double> &f, double uy0) const;

    /**
     * @param f vector containing distribution functions of a node
     * @param ux0 vertical velocity of the node
     * @return the pre-streaming vector of distribution functions of the node 
    */
    std::vector<double> calc_right_wall_bounce(std::vector<double> &f, double ux0) const;

    /**
     * @param f vector containing distribution functions of a node
     * @param ux0 vertical velocity of the node
     * @return the pre-streaming vector of distribution functions of the node 
    */
    std::vector<double> calc_bottom_wall_bounce(std::vector<double> &f, double ux0) const;


    std::vector<double> calc_top_left_corner_bounce() const;
    std::vector<double> calc_top_right_corner_bounce() const;
    std::vector<double> calc_bottom_left_corner_bounce() const;
    std::vector<double> calc_bottom_right_corner_bounce() const;

    virtual ~NEBB() = default;
};

#endif