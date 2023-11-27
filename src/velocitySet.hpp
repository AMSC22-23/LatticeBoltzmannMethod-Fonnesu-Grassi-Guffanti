#ifndef HH_VELOCITYSET
#define HH_VELOCITYSET


#include <iostream>
#include <string>
#include <vector>

struct WeightedDirection{
    std::vector<int> direction;
    double weight;
};

class VelocitySet
{
private:
    // TODO: vedere se cambiare implementazione

    std::size_t D;
    int Q;
    std::vector<WeightedDirection> Set;
    
public:
    VelocitySet()=default;
    
    const size_t get_dimensions() const{
        return D;
    }

    const int get_set_size() const{
        return Q;
    }

    const std::vector<WeightedDirection> get_velocity_set() const{
        return this->Set;
    }

    void set_dimension(const size_t D){
        this->D=D;
    }

    void set_q(const int Q){
        this->Q=Q;
    }

    void set_velocity_set();

    void initialize(const size_t D,const int Q){
        set_dimension(D);
        set_q(Q);
        set_velocity_set();
    }
    
};




#endif
