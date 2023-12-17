#include <iostream>
#include "lbm.hpp"

int main()
{
    lbm Prova(2, 1000.0, "../resources/lattices/lid_driven_cavity","TRT", 100);
    Prova.compute(20.0);

    return 0;
}
