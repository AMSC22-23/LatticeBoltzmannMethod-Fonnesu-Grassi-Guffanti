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
    
    // array of equilibrium populations (feq)
    std::vector <double> eq_populations;
    
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
    eq_populations ({0.0}),
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
     * Computes, returns and updates the macroscopic density scalar and velocity vector.
     * @param velocity_set the velocity set from which velocities are taken
     * @return the macroscopic density of the fluid and the velocity as a vector
    */
    const std::tuple<double, std::array<double, dim>> macroscopic_quantities(const VelocitySet& velocity_set) 
    {

        rho = 0.0;
        u = {0.0};
        // calculating the density
        for (const auto& population : populations) 
        {
            rho += population;
        }

        const std::vector<WeightedDirection> vs_content = velocity_set.get_velocity_set();
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

    /**
     * Initializes a lattice node by setting macroscopic density, velocity, populations and equilibrium populations
     * @param set_of_weights the weighted directions composing the velocity set
    */
    void initialize_node(const std::vector<WeightedDirection>& set_of_weights)
    {
        populations.resize(set_of_weights.size());
        eq_populations.resize(set_of_weights.size());

        if (type == FLUID)
        {
            // setting the macroscopic density
            rho = 1.0;

            // setting the macroscopic velocity
            u = {0.0};
            
            // setting the populations
            populations = {0.0};
            
            // setting the equilibrium populations: since u = 0 and rho = 1, the equation describing the 
            // equilibrium populations results only into the wi factor, which is the weight of a specific direction.
            std::size_t size = set_of_weights.size();
            for (std::size_t velocity_index = 0; velocity_index < size; ++velocity_index)
            {
                eq_populations[velocity_index] = set_of_weights[velocity_index].weight;
            }
        } else 
        {
            rho = 0.0;
            u = {0.0};
            populations = {0.0};
            eq_populations = {0.0};
        }
        
    }

    void compute_equilibrium_populations(const std::vector<WeightedDirection>& set_of_weights) 
    {
        //TODO: not yet implemented.
    }
};

#endif  // HH_LBM_LATTICE_NODE