#ifndef HH_LBM_LATTICE_NODE
#define HH_LBM_LATTICE_NODE

#include <array>
#include <tuple>
#include <algorithm>
#include "velocitySet.hpp"
#include "collisionModel.hpp"
#include "boundary.hpp"

/**
 * Enumeration that describes the type of node
*/
enum Node_type{FLUID, OPEN_BOUNDARY, BOUNDARY, SOLID}; 

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
    
    // array of populations after collision (fi*)
    std::vector<double> collision_populations;

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

    using PointInGrid = std::array<std::size_t, dim>;

    LatticeNode(
        const Node_type type_ = FLUID
    ) :
    type (type_),
    rho (0.0)
    {
        u.fill(0.0);
        collision_populations.resize(0);
        populations.resize(0);
        eq_populations.resize(0);
    };
    
    ~LatticeNode() = default;

    const std::vector<double> get_populations() const 
    {
        return populations;
    }

    const std::vector<double> get_eq_populations() const 
    {
        return eq_populations;
    }

    const std::vector<double> get_collision_populations() const
    {
        return collision_populations;
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

    const bool is_open_boundary() const 
    {
        return type == OPEN_BOUNDARY;
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

    double& set_population(std::size_t i)
    {
        return populations[i];
    }

    std::vector<double>& set_collision_populations()
    {
        return collision_populations;
    }

    std::array<double, dim>& set_u()
    {
        return u;
    }

    double& set_rho() 
    {
        return rho;
    }

    const std::tuple<double, std::array<double, dim>> get_macroscopic_quantities() const
    {
        return {rho, u};
    }

    /**
     * Computes, returns and updates the macroscopic density scalar and velocity vector.
     * @param velocity_set the velocity set from which velocities are taken
     * @return the macroscopic density of the fluid and the velocity as a vector
    */
    void update_macroscopic_quantities(const VelocitySet& velocity_set) 
    {
        rho = 0.0;
        u.fill(0.0);
        // calculating the density
        for (const auto& population : populations) 
        {
            rho += population;
        }

        const auto directions = velocity_set.get_velocity_set().direction;
        const std::size_t size = directions.size();

        // calculating the momentum density 
        for (std::size_t velocity_index = 0; velocity_index < size; ++velocity_index)
        {
            for (std::size_t i = 0; i < dim; ++i) 
            {
                u[i] += directions[velocity_index][i] * populations[velocity_index];
            }
        }

        // calculating the velocity by dividing the momentum by the density
        for (auto& momentum_density : u) {
            momentum_density = momentum_density / rho;
        }
    }

    /**
     * Initializes a lattice node by setting macroscopic density, velocity, populations and equilibrium populations
     * @param set_of_weights the weighted directions composing the velocity set
    */
    void initialize_fluid_node(const std::vector<double>& set_of_weights, const std::array<double, dim> node_u, const double node_rho = 1.0)
    {
        populations.resize(set_of_weights.size());
        eq_populations.resize(set_of_weights.size());
        collision_populations.resize(set_of_weights.size());

        // setting the macroscopic density
        rho = node_rho;

        // setting the macroscopic velocity
        u = node_u;
        
        // setting the populations
        std::fill(populations.begin(), populations.end(), 0.0);
        std::fill(collision_populations.begin(), collision_populations.end(), 0.0);
        
        // setting the equilibrium populations: since u = 0 and rho = 1, the equation describing the 
        // equilibrium populations results only into the wi factor, which is the weight of a specific direction.
        std::size_t size = set_of_weights.size();
        for (std::size_t velocity_index = 0; velocity_index < size; ++velocity_index)
        {
            eq_populations[velocity_index] = set_of_weights[velocity_index];
        }
    }

    /**
     * Initializes an open boundary node by setting its parameters based on functions passed as input
     * @param velocity_set the velocity set, used to extract weights.
     * @param rho_w desnity of the node
     * @param u_w velocities of the node
     * @param one_over_speed_of_sound_squared the reciprocal of the square of the speed of sound of the system. Defaults to 3.0.
    */
    void initialize_open_boundary_node(const VelocitySet& velocity_set, const double rho_w, const std::array<double, dim>& u_w, const double one_over_speed_of_sound_squared = 3.0)
    {
        rho = rho_w;
        u = u_w;
        const auto directions = velocity_set.get_velocity_set().direction;
        std::vector<double> weights = velocity_set.get_velocity_set().weight;
        
        collision_populations.resize(weights.size());
        eq_populations.resize(weights.size());
        populations.resize(weights.size());

        std::fill(eq_populations.begin(), eq_populations.end(), 0.0);
        std::fill(collision_populations.begin(), collision_populations.end(), 0.0);

        const std::size_t size = directions.size();

        for(std::size_t i=0; i < size; ++i)
        {
            double cu = directions[i][0] * u[0] + directions[i][1] * u[1];
            populations[i] = weights[i] * rho * (
                1.0 + one_over_speed_of_sound_squared * cu + 
                0.5 * one_over_speed_of_sound_squared * one_over_speed_of_sound_squared * cu * cu - 
                0.5 * one_over_speed_of_sound_squared * (u[0] * u[0] + u[1] * u[1]));
        }

        return;
    }

    /**
     * Initializes a generic boundary by setting all the parameters to 0.0
    */
    void initialize_generic_boundary(const VelocitySet& velocity_set)
    {   
        
        rho = 0.0;
        u.fill(0.0);
        
        std::size_t size = velocity_set.get_set_size();

        collision_populations.resize(size);
        eq_populations.resize(size);
        populations.resize(size);

        std::fill(eq_populations.begin(), eq_populations.end(), 0.0);
        std::fill(populations.begin(), populations.end(), 0.0);
        std::fill(collision_populations.begin(), collision_populations.end(), 0.0);

    }

    /**
     * Method that computes the equilibrium populations using 
     * a given speed of sound and a set. The speed of sound may be different due
     * to debugging purposes or as DeltaT/Tau is not 1
     * @param velocity_set the velocity set, used to extract weights.
     * @param one_over_speed_of_sound_squared the reciprocal of the square of the speed of sound of the system. Defaults to 3.0.
    */
    void compute_equilibrium_populations(const WeightedDirection& velocity_set, const double one_over_speed_of_sound_squared = 3.0) 
    {   
        double u_dot_u;
        double u_dot_ci;

        auto directions = velocity_set.direction;
        auto weights = velocity_set.weight;
        auto size = directions.size();

        for (std::size_t i = 0; i < size; ++i)
        {
            // computing the dot products
            u_dot_u = 0.0;
            u_dot_ci = 0.0;
            for (std::size_t d = 0; d < dim; d++)
            {
                u_dot_u +=  u[d] * u[d];
                u_dot_ci += u[d] * directions[i][d];
            }

            eq_populations[i] = weights[i] * rho * (
                1.0 + one_over_speed_of_sound_squared * u_dot_ci + 
                0.5 * one_over_speed_of_sound_squared * one_over_speed_of_sound_squared * u_dot_ci * u_dot_ci -
                0.5 * one_over_speed_of_sound_squared * u_dot_u);
        }
    }
};

#endif  // HH_LBM_LATTICE_NODE