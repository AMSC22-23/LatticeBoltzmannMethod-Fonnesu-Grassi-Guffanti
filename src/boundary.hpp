#ifndef HH_BOUNDARY
#define HH_BOUNDARY

#include <vector>

class Boundary
{
public: 
    virtual std::vector<double> calcBounce(const std::vector<double> &f, const std::vector<double> &fEq, double tau, double deltaT) const = 0;
    virtual ~Boundary() = default;
};

class BounceBack : Boundary
{
public:
    virtual std::vector<double> calcBounce(const std::vector<double> &f, const std::vector<double> &fEq, double tau, double deltaT) const override;
};

#endif