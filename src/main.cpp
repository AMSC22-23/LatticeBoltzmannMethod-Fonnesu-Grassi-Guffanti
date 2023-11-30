#include <iostream>
#include "lbm.hpp"

int main()
{
    lbm Prova(2,9,"../resources/lattices/lid_driven_cavity/2d_100_100_lid_driven","BGK","ZH");
    Prova.print();
    /*
    PROVA CREZIONE VELOCITYSETS
    VelocitySet ProvaD1Q3(1,3);
    VelocitySet ProvaD2Q9(2,9);
    VelocitySet ProvaD3Q15(3,15);
    VelocitySet ProvaD3Q19(3,19);
    VelocitySet ProvaD3Q27(3,27);*/
    return 0;
}
