#ifndef HH_INLET_INITIALIZER
#define HH_INLET_INITIALIZER

#include <array>
#include <cmath>
#include <iostream>

/**
 * Abstract class that represents the object capable of initializing boundary nodes in a lattice
 * by setting the values of the inlet fields. The object is generic with respect to the number of 
 * dimensions. We provide some default lattice initializers. 
 * @author Luca Guffanti
*/
template<std::size_t dim>
class InletInitializer
{
public:
    InletInitializer() = default;
    virtual ~InletInitializer() = default;

    
    /**
     * Calculates the velocity field and density field of a boundary node 
     * knowing its position in the lattice and the iteration of the simulation
     * @param point coordinates of a point in the grid
     * @param it number of the iteration
     * @return a tuple containing the new velocity of the point
    */
    virtual std::tuple<std::array<double, dim>, double> set_inlets(
        const std::array<std::size_t, dim>& point,
        const std::size_t it) const = 0;
};

/**
 * This initializer sets the fields of the top wall in the lid driven cavity simulation.
 * @author Luca Guffanti
*/
class LidDrivenCavityUniformInitializer : public InletInitializer<2>
{
private:
    /**
     * Constant used to calculate the fields
    */
    const double sigma;
    /**
     * Scaling factor for the calculated vertical velocity
    */
    const double u_y_scaling;

    /**
     * Width of the computational domain
    */
    const double width;

public:
    LidDrivenCavityUniformInitializer(
        const std::size_t width_,
        const double u_y_scaling_
    ) :
    sigma (10 * width_),
    u_y_scaling (u_y_scaling_),
    width (width_)
    {};
    
    virtual ~LidDrivenCavityUniformInitializer() override = default;

    virtual std::tuple<std::array<double, 2>, double> set_inlets(
        const std::array<std::size_t, 2>& point,
        const std::size_t it) const override
    {
        double ret  = (1.0 - std::exp(-(it * it * 1.0)/(2.0*sigma*sigma)));
        if (point[0] == 0)
        {
            //std::cout << ret << std::endl;
            return {{ret * 0.2, 0.0}, 0.0};
        } else
        {
            return {{0.0, 0.0}, 0.0};
        }
    }
};


#endif // HH_INLET_INITIALIZER