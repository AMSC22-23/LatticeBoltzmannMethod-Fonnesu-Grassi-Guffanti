#ifndef HH_VELOCITYSET
#define HH_VELOCITYSET


#include <iostream>
#include <string>
#include <vector>


class VelocitySet
{
private:
    // TODO: vedere se cambiare implementazione

    const std::size_t D;
    const int Q;
    std::vector<std::vector<int>> directions;
    std::vector<double> weigths;
public:
    VelocitySet(const std::size_t D, const int Q);
    

    int getDirections(int i);
    int getWeights(int i);
};


#endif
