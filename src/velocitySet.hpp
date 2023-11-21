#ifndef HH_VELOCITYSET
#define HH_VELOCITYSET


#include <iostream>
#include <string>
#include <vector>

class VelocitySet
{
private:
    // TODO: vedere se cambiare implementazione

    int D;
    int Q;
    std::vector<int> directions;
    std::vector<int> weigths;
public:
    VelocitySet(std::string velocitySet);

    int getDirection(int i);
    int getWeight(int i);
};

#endif
