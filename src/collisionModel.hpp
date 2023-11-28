#ifndef HH_COLLISION_MODEL
#define HH_COLLISION_MODEL

#include <vector>

/**
 * General class that will be extended by other classes to implement vatious collision models
 * @author Lorenzo Fonnesu
*/
class CollisionModel
{
public:
    virtual std::vector<double> calc_collision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const = 0;

    virtual ~CollisionModel() = default;
};

/**
 * Implementation of the BGK method for collisions
 * @author Lorenzo Fonnesu
*/
class BGK : public CollisionModel
{
public:
    /**
    * @param f vector containing distribution functions of a node
    * @param fEq vector containing equilibrium distribution functions of a node
    * @param tConst deltaT/tau
    * @param tConj 1-tConst
    * @return the post-collision distribution functions
    */
    virtual std::vector<double> calc_collision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;

    virtual ~BGK() = default; 
};

class TRT : public CollisionModel
{
public:
    virtual std::vector<double> calc_collision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;

    virtual ~TRT() = default;
};

class MRT : public CollisionModel
{
public:
    std::vector<double> calc_collision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;

    virtual ~MRT() = default;
};

#endif