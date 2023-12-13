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
    BGK() = default;
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
private: 
    std::vector<int> f_conj = {0, 3, 4, 1, 2, 7, 8, 5, 6}; //Vector containing the opposite index for each direction (Example: 1<->3, 5<->7)
    const double om_p = 1; //1 / tau
    const double om_m = 1; //1 / ((0.25 / 1 - 0.5) + 0.5)
public:
    TRT() = default;
    virtual std::vector<double> calc_collision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;

    virtual ~TRT() = default;
};

class MRT : public CollisionModel
{
public:
    MRT() = default;
    std::vector<double> calc_collision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;

    virtual ~MRT() = default;
};

#endif