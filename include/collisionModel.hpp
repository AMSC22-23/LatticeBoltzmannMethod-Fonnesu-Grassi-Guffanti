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
    /**
    * @param f vector containing distribution functions of a node
    * @param fEq vector containing equilibrium distribution functions of a node
    * @return the post-collision distribution functions
    */
    virtual std::vector<double> calc_collision(const std::vector<double> &f, const std::vector<double> &fEq) const = 0;

    virtual ~CollisionModel() = default;
};

/**
 * Implementation of the BGK method for collisions
 * @author Lorenzo Fonnesu
*/
class BGK : public CollisionModel
{

private:
    const double t_conj;
    const double t_const;

public:
    /** 
    * @param tConst deltaT/tau
    * @param tConj 1-tConst
    */
    BGK(const double t_const_, const double t_conj_);

    virtual std::vector<double> calc_collision(const std::vector<double> &f, const std::vector<double> &fEq) const override;

    virtual ~BGK() = default; 
};

/**
 * Implementation of the TRT method for collision model
*/
class TRT : public CollisionModel
{
private: 
    double om_p; 
    double om_m;
    const double lambda = 1.0/4.0;
    const double tau_p;
    const double tau_m;
    //Vector containing the opposite index for each direction (Example: 1<->3, 5<->7)
    std::vector<int> f_conj = {0, 3, 4, 1, 2, 7, 8, 5, 6}; 

    public:
    /**
     * @param tau time constant
    */
    TRT(const double tau);

    virtual std::vector<double> calc_collision(const std::vector<double> &f, const std::vector<double> &fEq) const override;

    virtual ~TRT() = default;
};

class MRT : public CollisionModel
{
public:
    MRT() = default;
    std::vector<double> calc_collision(const std::vector<double> &f, const std::vector<double> &fEq) const override;

    virtual ~MRT() = default;
};

#endif