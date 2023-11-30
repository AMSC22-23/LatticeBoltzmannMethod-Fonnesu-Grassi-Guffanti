#include "boundary.hpp"

std::vector<double> NEBB::calc_top_wall_bounce(std::vector<double> &f, double uy0) const
{
    std::vector<double> fTemp = f;
    static constexpr double two_thirds = 2.0/3.0;
    static constexpr double one_sixth = 1.0/6.0;
    static constexpr double one_half = 0.5;

    double rho = (f[0] + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[6])) / (1.0 + uy0);

    double ru = rho * uy0;

    fTemp[4] = f[2] - two_thirds * ru;
    fTemp[7] = f[5] - one_sixth * ru + one_half * (f[1] - f[3]);
    fTemp[8] = f[6] - one_sixth * ru + one_half * (f[3] - f[1]);

    return fTemp;
}

std::vector<double> NEBB::calc_bottom_wall_bounce(std::vector<double> &f, double uy0) const
{
    std::vector<double> fTemp = f;
    static constexpr double two_thirds = 2.0/3.0;
    static constexpr double one_sixth = 1.0/6.0;
    static constexpr double one_half = 0.5;

    double rho = (f[0] + f[1] + f[3] + 2.0 * (f[4] + f[7] + f[8])) / (1.0 - uy0);

    double ru = rho * uy0;

    fTemp[2] = f[4] + two_thirds * ru;
    fTemp[5] = f[7] + one_sixth * ru - one_half * (f[1] - f[3]);
    fTemp[6] = f[8] + one_sixth * ru - one_half * (f[3] - f[1]);

    return fTemp;
}

std::vector<double> NEBB::calc_left_wall_bounce(std::vector<double> &f, double ux0) const
{
    std::vector<double> fTemp = f;
    static constexpr double two_thirds = 2.0/3.0;
    static constexpr double one_sixth = 1.0/6.0;
    static constexpr double one_half = 0.5;

    double rho = (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[7] + f[6])) / (1.0 - ux0);

    double ru = rho * ux0;

    fTemp[1] = f[3] + two_thirds * ru;
    fTemp[5] = f[7] + one_sixth * ru - one_half * (f[2] - f[4]);
    fTemp[8] = f[6] + one_sixth * ru - one_half * (f[4] - f[2]);

    return fTemp;
}

std::vector<double> NEBB::calc_right_wall_bounce(std::vector<double> &f, double ux0) const
{
    std::vector<double> fTemp = f;
    static constexpr double two_thirds = 2.0/3.0;
    static constexpr double one_sixth = 1.0/6.0;
    static constexpr double one_half = 0.5;

    double rho = (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / (1.0 + ux0);

    double ru = rho * ux0;

    fTemp[3] = f[1] - two_thirds * ru;
    fTemp[7] = f[5] - one_sixth * ru + one_half * (f[2] - f[4]);
    fTemp[6] = f[8] - one_sixth * ru + one_half * (f[4] - f[2]);

    return fTemp;

}

std::vector<double> NEBB::calc_top_left_corner_bounce() const
{

}

std::vector<double> NEBB::calc_top_right_corner_bounce() const
{

}

std::vector<double> NEBB::calc_bottom_left_corner_bounce() const
{

}

std::vector<double> NEBB::calc_bottom_right_corner_bounce() const
{

}