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
    /**
    * @param f array containing distribution functions of a node
    * @param fEq array containing equilibrium distribution functions of a node
    * @param tConst deltaT/tau
    * @param tConj 1-tConst
    * @return the post-bounce distribution functions
    */
    std::vector<double> calcBounce(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const override;

    /*
    std::array<double, q> calcBounce(const std::array<double, q> &f, const std::array<double, q> &fEq, double tConst, double tConj) const override
    {
        std::array<double, q> fStar;

        for(size_t i=0; i<q; i++){
            // The idea is to look, with the help of c values and a getter of single nodes in the lattice (maybe the coord of the current
            // node too), if the current node neighbours are fluid or not
            // If they are, we can proceed as BGK, otherwise following the rules of BB method
            // Another option is to call this method singularly for the nodes that aren't fluid
            // Finally, please note that there are some cases (for example corners) in which this method doesn't work in the same way
        }
    }
    */

    virtual ~BounceBack() = default;
};

#endif