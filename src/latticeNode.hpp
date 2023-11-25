#ifndef HH_LBM_LATTICE_NODE
#define HH_LBM_LATTICE_NODE

#include <array>
#include <tuple>
#include <algorithm>
#include "velocitySet.hpp"

/**
 * Enumeration that describes the type of node
*/
enum Node_type{FLUID, BOUNDARY, SOLID};

/**
 * object that represents a single node of the lattice:
 * each lattice node can either be a solid, boundary or fluid node, and it
 * contains information regarding macroscopic quantities and the populations.
 * 
 * @note templating is used to generalize the number of dimensions and velocities
 * 
 * @author Luca Guffanti
*/
template <std::size_t dim>
class LatticeNode
{
private:
    // array of discrete populations (fi)
    std::vector<double> populations;
    
    // type of the node
    Node_type type;
    
    // macroscopic quantities
    // density of the fluid
    double rho;
    // velocity of the fluid
    std::array<double, dim> u;


public:
    LatticeNode(
        const Node_type type_ = FLUID
    ) :
    type (type_),
    populations ({0.0}),
    rho (0.0),
    u ({0.0})
    {};
    
    ~LatticeNode() = default;

    const std::vector<double> get_populations() const 
    {
        return populations;
    }

    const std::array<double, dim> get_u() const 
    {
        return u;
    }

    const double get_rho() const 
    {
        return rho;
    }

    const bool is_fluid() const 
    {
        return type == FLUID;
    }

    const bool is_boundary() const 
    {
        return type == BOUNDARY;
    }
    
    const bool is_solid() const 
    {
        return type == SOLID;
    }

    Node_type& set_type()
    {
        return type;
    }

    std::vector<double>& set_populations()
    {
        return populations;
    }

    std::array<double, dim>& set_u()
    {
        return u;
    }

    double& set_rho() 
    {
        return rho;
    }


    /**
     * @param velocity_set the velocity set from which velocities are taken
     * @return the macroscopic density of the fluid and the velocity as a vector
    */
    const std::tuple<double, std::array<double, dim>> macroscopic_quantities(const VelocitySet& velocity_set) 
    {

        rho = 0.0;
        u = {0.0};
        populations = {1.0, 2.0, 3.0};
        // calculating the density
        for (const auto& population : populations) 
        {
            rho += population;
        }

        const std::vector<WeightedDirection> vs_content = velocity_set.getVelocitySets();
        const std::size_t size = vs_content.size();

        // calculating the momentum density 
        for (std::size_t velocity_index = 0; velocity_index < size; ++velocity_index)
        {
            for (std::size_t i = 0; i < dim; ++i) 
            {
                u[i] += vs_content[velocity_index].direction[i] * populations[velocity_index];
            }
        }

        // calculating the velocity by dividing the momentum by the density
        for (auto& momentum_density : u) {
            momentum_density = momentum_density / rho;
        }
        return {rho, u};
    }

    void initialize_node(const double rho, const std::array<double, dim> u)
    {
        
    }
};

#endif  // HH_LBM_LATTICE_NODE