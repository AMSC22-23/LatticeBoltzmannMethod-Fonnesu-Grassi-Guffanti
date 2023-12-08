#include "boundary.hpp"

void NEBB::calc_top_wall_bounce(LatticeNode<2> &node) const
{ 
    std::vector<double> f = node.get_populations();

    const double rho = (f[0] + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[7])) / (1.0 + node.get_u()[1]);

    node.set_population(4) = f[2] - two_thirds*rho*node.get_u()[1];

    node.set_population(8) = f[6] - one_half*(f[1] - f[3]) +
                             one_half*rho*node.get_u()[0] -
                             one_sixth*rho*node.get_u()[1];

    node.set_population(7) = f[5] + one_half*(f[1] - f[3]) -
                             one_half*rho*node.get_u()[0] -
                             one_sixth*rho*node.get_u()[1];
}

void NEBB::calc_bottom_wall_bounce(LatticeNode<2> &node) const
{
    std::vector<double> f = node.get_populations();

    const double rho = (f[0] + f[1] + f[3] + 2.0 * (f[4] + f[7] + f[8])) / (1.0 - node.get_u()[1]);

    const double ru = rho * node.get_u()[1];

    node.set_population(2) = f[4] + two_thirds*rho*node.get_u()[1];

    node.set_population(5) = f[7] - one_half*(f[1]-f[3]) + 
                             one_half*rho*node.get_u()[0] + 
                             one_sixth*rho*node.get_u()[1];

    node.set_population(6) = f[8] + one_half*(f[1]-f[3]) - 
                             one_half*rho*node.get_u()[0] + 
                             one_sixth*rho*node.get_u()[1];
}

void NEBB::calc_left_wall_bounce(LatticeNode<2> &node) const
{
    std::vector<double> f = node.get_populations();

    const double rho = (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[7] + f[6])) / (1.0 - node.get_u()[0]);

    const double ru = rho * node.get_u()[0];

    node.set_population(1) = f[3] + two_thirds*rho*node.get_u()[0];

    node.set_population(5) = f[7] - one_half*(f[2]-f[4]) + 
                             one_sixth*rho*node.get_u()[0] + 
                             one_half*rho*node.get_u()[1];

    node.set_population(8) = f[6] + one_half*(f[2]-f[4]) + 
                             one_sixth*rho*node.get_u()[0] - 
                             one_half*rho*node.get_u()[1];
}

void NEBB::calc_right_wall_bounce(LatticeNode<2> &node) const
{
    std::vector<double> f = node.get_populations();

    const double rho = (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / (1.0 + node.get_u()[0]);

    const double ru = rho * node.get_u()[0];

    node.set_population(3) = f[1] - two_thirds*rho*node.get_u()[0];

    node.set_population(7) = f[5] + one_half*(f[2]-f[4]) - 
                             one_sixth*rho*node.get_u()[0] - 
                             one_half*rho*node.get_u()[1];

    node.set_population(6) = f[8] - one_half*(f[2]-f[4]) - 
                             one_sixth*rho*node.get_u()[0] + 
                             one_half*rho*node.get_u()[1];
}

void NEBB::calc_top_left_corner_bounce(LatticeNode<2> &node) const
{
    node.set_population(1) = node.get_populations()[3] + two_thirds * node.get_rho() * node.get_u()[0];
    node.set_population(4) = node.get_populations()[2] - two_thirds * node.get_rho() * node.get_u()[1];
    node.set_population(8) = node.get_populations()[6] + one_sixth * node.get_rho() * (node.get_u()[0] - node.get_u()[1]);

    node.set_population(5) = 0.0;
    node.set_population(7) = 0.0;

    node.set_population(0) = node.get_rho() 
                            - node.get_populations()[1] - node.get_populations()[3] - node.get_populations()[2] - node.get_populations()[4]
                            - node.get_populations()[5] - node.get_populations()[7] - node.get_populations()[6] - node.get_populations()[8]; 
}

void NEBB::calc_top_right_corner_bounce(LatticeNode<2> &node) const
{
    node.set_population(3) = node.get_populations()[1] - two_thirds * node.get_rho() * node.get_u()[0];
    node.set_population(4) = node.get_populations()[3] - two_thirds * node.get_rho() * node.get_u()[1];
    node.set_population(7) = node.get_populations()[5] - one_sixth * node.get_rho() * (node.get_u()[0] + node.get_u()[1]);

    node.set_population(6) = 0.0;
    node.set_population(8) = 0.0;

    node.set_population(0) = node.get_rho() 
                            - node.get_populations()[1] - node.get_populations()[3] - node.get_populations()[2] - node.get_populations()[4]
                            - node.get_populations()[5] - node.get_populations()[7] - node.get_populations()[6] - node.get_populations()[8];
}

void NEBB::calc_bottom_left_corner_bounce(LatticeNode<2> &node) const
{
    node.set_population(1) = node.get_populations()[3] - two_thirds * node.get_rho() * node.get_u()[0];
    node.set_population(2) = node.get_populations()[4] - two_thirds * node.get_rho() * node.get_u()[1];
    node.set_population(5) = node.get_populations()[7] + one_sixth * node.get_rho() * (node.get_u()[0] + node.get_u()[1]);

    node.set_population(6) = 0.0;
    node.set_population(8) = 0.0;

    node.set_population(0) = node.get_rho() 
                            - node.get_populations()[1] - node.get_populations()[2] - node.get_populations()[3] - node.get_populations()[4]
                            - node.get_populations()[5] - node.get_populations()[6] - node.get_populations()[7] - node.get_populations()[8];
}

void NEBB::calc_bottom_right_corner_bounce(LatticeNode<2> &node) const
{
    node.set_population(3) = node.get_populations()[1] - two_thirds * node.get_rho() * node.get_u()[0];
    node.set_population(2) = node.get_populations()[4] + two_thirds * node.get_rho() * node.get_u()[1];
    node.set_population(6) = node.get_populations()[8] - one_sixth * node.get_rho() * (node.get_u()[0] + node.get_u()[1]);

    node.set_population(5) = 0.0;
    node.set_population(7) = 0.0;

    node.set_population(0) = node.get_rho() 
                            - node.get_populations()[1] - node.get_populations()[2] - node.get_populations()[3] - node.get_populations()[4]
                            - node.get_populations()[5] - node.get_populations()[6] - node.get_populations()[7] - node.get_populations()[8];
}