#include "collisionModel.hpp"

// ===================BGK===================
BGK::BGK(const double t_conj_, const double t_const_):
t_conj (t_conj),
t_const (t_const)
{}

std::vector<double> BGK::calc_collision(const std::vector<double> &f, const std::vector<double> &fEq) const
{
    std::vector<double> fStar;
    fStar.resize(f.size());
    // This loop works for all DxQy, we can improve it deciding to remove D3Q19, because all the others y can be divided by 3
    for(std::size_t i=0; i<f.size(); i++){
        fStar[i] = f[i] * t_conj + fEq[i] * t_const;
    }

    return fStar;
}

// ===================TRT===================
TRT::TRT(const double tau):
tau_p (tau),
tau_m (lambda / (tau_p - 0.5) + 0.5)
{
    om_p = 1/tau_p;
    om_m = 1/tau_m;
}

//Others models, we can decide to use them later
std::vector<double> TRT::calc_collision(const std::vector<double> &f, const std::vector<double> &fEq) const
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

// ===================MRT===================
std::vector<double> MRT::calc_collision(const std::vector<double> &f, const std::vector<double> &fEq) const
{
    return {0};
}
