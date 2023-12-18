#include <iostream>
#include "lbm.hpp"

int main()
{
    lbm Prova(2,100.0, "../resources/lattices/lid_driven_cavity/weak_scalability/100","BGK", 100);
    //Prova.compute(20.0);
    Prova.perform_weak_scaling_test(100);
    return 0;
}
