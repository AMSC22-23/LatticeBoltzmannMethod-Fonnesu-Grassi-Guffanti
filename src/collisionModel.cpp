#include "collisionModel.hpp"

std::vector<double> BGK::calc_collision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const
{
    std::vector<double> fStar;

    // This loop works for all DxQy, we can improve it deciding to remove D3Q19, because all the others y can be divided by 3
    for(int i=0; i<f.size(); i++){
        fStar[i] = f[i] * tConj + fEq[i] * tConst;
    }

    return fStar;
}

//Others models, we can decide to use them later
std::vector<double> TRT::calc_collision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const
{

}

std::vector<double> MRT::calc_collision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const
{

}
