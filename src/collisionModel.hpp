#ifndef HH_COLLISION_MODEL
#define HH_COLLISION_MODEL

#include <vector>
//#include <array>

//template <std::size_t q>
class CollisionModel
{
public:
    virtual std::vector<double> calcCollision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const = 0;

    //virtual std::array<double, q> calcCollision(const std::array<double, q> &f, const std::array<double, q> &fEq, double tConst, double tConj) const = 0;

    virtual ~CollisionModel() = default;
};

//template <std::size_t q>
class BGK : public CollisionModel //<q>
{
public:
    /**
    * @param f array containing distribution functions of a node
    * @param fEq array containing equilibrium distribution functions of a node
    * @param tConst deltaT/tau
    * @param tConj 1-tConst
    * @return the post-collision distribution functions
    */
    std::vector<double> calcCollision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;

    /*
    std::array<double, q> calcCollision(const std::array<double, q> &f, const std::array<double, q> &fEq, double tConst, double tConj) const override
    {
        std::array<double, q> fStar;

        // This loop works for all DxQy, we can improve it deciding to remove D3Q19, because all the others y can be divided by 3
        for(size_t i=0; i<q; i++){
            fStar[i] = f[i] * tConj + fEq[i] * tConst;
        }

        return fStar;
        // To assign correctly all the new distribution functions, we can use c values
    }
    */

    virtual ~BGK() = default; 
};

class TRT : public CollisionModel
{
public:
    std::vector<double> calcCollision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;

    /*
    std::array<double, q> calcCollision(const std::array<double, q> &f, const std::array<double, q> &fEq, double tConst, double tConj) const override
    {
        
    }
    */

    virtual ~TRT() = default;
};

class MRT : public CollisionModel
{
public:
    std::vector<double> calcCollision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;

    /*
    std::array<double, q> calcCollision(const std::array<double, q> &f, const std::array<double, q> &fEq, double tConst, double tConj) const override
    {
        
    }
    */

    virtual ~MRT() = default;
};

#endif