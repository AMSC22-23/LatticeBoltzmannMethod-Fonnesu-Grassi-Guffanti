#ifndef HH_BOUNDARY
#define HH_BOUNDARY

#include <vector>
//#include <array>

//template <std::size_t q>
class Boundary
{
public: 
    virtual std::vector<double> calcBounce(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const = 0;

    //virtual std::array<double, q> calcBounce(const std::array<double, q> &f, const std::array<double, q> &fEq, double tConst, double tConj) const = 0;

    virtual ~Boundary() = default;
};

//template <std::size_t q>
class BounceBack : Boundary //<q>
{
public:
    /*
    * @param f array containing distribution functions of a node
    * @param fEq array containing equilibrium distribution functions of a node
    * @param tConst deltaT/tau
    * @param tConj: 1-tConst
    * @return the post-bounce distribution functions
    */
    std::vector<double> calcBounce(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;

    /*
    std::array<double, q> calcBounce(const std::array<double, q> &f, const std::array<double, q> &fEq, double tConst, double tConj) const override
    {
        
    }
    */

    virtual ~BounceBack() = default;
};

#endif