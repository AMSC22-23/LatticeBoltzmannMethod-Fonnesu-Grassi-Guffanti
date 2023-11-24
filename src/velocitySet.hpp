#ifndef HH_VELOCITYSET
#define HH_VELOCITYSET


#include <iostream>
#include <string>
#include <vector>


class VelocitySet
{
private:
    // TODO: vedere se cambiare implementazione

    struct WeightedDirection{
        std::vector<int> direction;
        double weight;
    };

    const std::size_t D;
    const int Q;
    std::vector<WeightedDirection> Set;
    
public:
    VelocitySet(const std::size_t D, const int Q);
    
    const size_t getDimensions() const{
        return D;
    }

    const int getSetSize() const{
        return Q;
    }

    const std::vector<WeightedDirection> getVelocitySets() const{
        return this->Set;
    }
};




#endif
