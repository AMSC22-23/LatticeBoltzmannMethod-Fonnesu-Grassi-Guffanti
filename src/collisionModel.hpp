#ifndef HH_COLLISION_MODEL
#define HH_COLLISION_MODEL

#include <vector>

class CollisionModel
{
public:
    // Virtual function to calculate f^*. It will be implemented by all the collision models
    virtual std::vector<double> calcCollision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const = 0;
    virtual ~CollisionModel() = default;
};

class BGK : public CollisionModel
{
public:
    std::vector<double> calcCollision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;
};

class TRT : public CollisionModel
{
public:
    std::vector<double> calcCollision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;
};

class MRT : public CollisionModel
{
public:
    std::vector<double> calcCollision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;
};

#endif