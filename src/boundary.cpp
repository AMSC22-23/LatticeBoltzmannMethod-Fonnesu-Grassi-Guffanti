#include "boundary.hpp"

std::vector<double> NEBB::calc_top_wall_bounce(std::vector<double> &f, double uy0) const
{
    std::vector<double> fTemp = f;

    double rho = (f[0] + f[1] + f[3] + 2 * (f[2] + f[5] + f[6])) / (1 + uy0);

    double ru = rho * uy0;

    fTemp[4] = f[2] - 2/3 * ru;
    fTemp[7] = f[5] - 1/6 * ru + 1/2 * (f[1] - f[3]);
    fTemp[8] = f[6] - 1/6 * ru + 1/2 * (f[3] - f[1]);

    return fTemp;
}

std::vector<double> NEBB::calc_bottom_wall_bounce(std::vector<double> &f, double uy0) const
{
    std::vector<double> fTemp = f;

    double rho = (f[0] + f[1] + f[3] + 2 * (f[4] + f[7] + f[8])) / (1 - uy0);

    double ru = rho * uy0;

    fTemp[2] = f[4] + 2/3 * ru;
    fTemp[5] = f[7] + 1/6 * ru - 1/2 * (f[1] - f[3]);
    fTemp[6] = f[8] + 1/6 * ru - 1/2 * (f[3] - f[1]);

    return fTemp;
}

std::vector<double> NEBB::calc_left_wall_bounce(std::vector<double> &f, double ux0) const
{
    std::vector<double> fTemp = f;

    double rho = (f[0] + f[2] + f[4] + 2 * (f[3] + f[7] + f[6])) / (1 - ux0);

    double ru = rho * ux0;

    fTemp[1] = f[3] + 2/3 * ru;
    fTemp[5] = f[7] + 1/6 * ru - 1/2 * (f[2] - f[4]);
    fTemp[8] = f[6] + 1/6 * ru - 1/2 * (f[4] - f[2]);

    return fTemp;
}

std::vector<double> NEBB::calc_right_wall_bounce(std::vector<double> &f, double ux0) const
{
    std::vector<double> fTemp = f;

    double rho = (f[0] + f[2] + f[4] + 2 * (f[1] + f[5] + f[8])) / (1 + ux0);

    double ru = rho * ux0;

    fTemp[3] = f[1] - 2/3 * ru;
    fTemp[7] = f[5] - 1/6 * ru + 1/2 * (f[2] - f[4]);
    fTemp[6] = f[8] - 1/6 * ru + 1/2 * (f[4] - f[2]);

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