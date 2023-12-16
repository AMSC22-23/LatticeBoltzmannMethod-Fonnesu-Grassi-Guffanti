#include "collisionModel.hpp"

std::vector<double> BGK::calc_collision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const
{
    std::vector<double> fStar;
    fStar.resize(f.size());
    // This loop works for all DxQy, we can improve it deciding to remove D3Q19, because all the others y can be divided by 3
    for(std::size_t i=0; i<f.size(); i++){
        fStar[i] = f[i] * tConj + fEq[i] * tConst;
    }

    return fStar;
}

//Others models, we can decide to use them later
std::vector<double> TRT::calc_collision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const
{
    int j;
    std::vector<double> fStar;
    fStar.resize(f.size());

    fStar[0] = (1 - om_p) * f[0] + om_p * fEq[0];

    for(std::size_t i=1; i<f.size(); i++){
        j = f_conj[i];
        fStar[i] = (1 - 0.5 * (om_p + om_m)) * f[i] - 
                       0.5 * (om_p - om_m) * f[j] +
                       0.5 * (om_p + om_m) * fEq[i] +
                       0.5 * (om_p - om_m) * fEq[j];
    }

    return fStar;
}

TRT::TRT(const double tau):
    tau_p (tau),
    tau_m (lambda / (tau_p - 0.5) + 0.5)
    {
        om_p = 1/tau_p;
        om_m = 1/tau_m;
    }

std::vector<double> MRT::calc_collision(const std::vector<double> &f, const std::vector<double> &fEq, double tConst, double tConj) const
{
    return {0};
}
