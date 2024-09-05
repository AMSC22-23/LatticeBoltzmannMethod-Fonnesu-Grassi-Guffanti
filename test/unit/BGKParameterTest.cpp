#include "llalbm.hpp"
#include <iostream>

int main()
{
    using namespace llalbm::core;

    collisions::BGKCollisionPolicy<2>::initialize(0.51, 0.01);
    collisions::BGKCollisionPolicy<2>::allowed_tau(0.02, 10);
    return 0;
}