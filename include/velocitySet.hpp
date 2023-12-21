#ifndef HH_VELOCITYSET
#define HH_VELOCITYSET


#include <iostream>
#include <string>
#include <vector>

/**
 * The velocity set is described here as a structure of array. The direction vector contains 
 * the direction of the velocities pointing towards neighbouring nodes, while the weight vector contains the
 * respective weights.
*/
struct WeightedDirection{
    std::vector<std::vector<int>> direction;
    std::vector<double> weight;
};

/**
 * The data structure representing a velocity set: it contains the number of dimensions, the number
 * of directions, and the set of directions coupled with weights. 
*/
class VelocitySet
{
private:

    std::size_t D;
    std::size_t Q;
    WeightedDirection Set;
    
public:
    VelocitySet()=default;
    
    const size_t get_dimensions() const{
        return D;
    }

    const std::size_t get_set_size() const{
        return Q;
    }

    const WeightedDirection get_velocity_set() const{
        return this->Set;
    }

    void set_dimension(const size_t D){
        this->D=D;
    }

    void set_q(const int Q){
        this->Q=Q;
    }

    /**
     * Reads the velocity set file and populates the velocity set
    */
    void set_velocity_set();

    /**
     * Initializes the velocity set
    */
    void initialize(const std::size_t D, const std::size_t Q){
        set_dimension(D);
        set_q(Q);
        set_velocity_set();
    }
    
};




#endif
