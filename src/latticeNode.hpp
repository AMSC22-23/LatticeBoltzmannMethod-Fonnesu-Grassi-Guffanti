#ifndef HH_LBM_LATTICE_NODE
#define HH_LBM_LATTICE_NODE

#include <array>
#include <tuple>

#include "velocitySet.hpp"

template <std::size_t dim, std::size_t numVelocities>
class LatticeNode
{
private:
    // array of discrete populations (fi)
    std::array<double, numVelocities> populations;
    
    // whether a lattice node is a fluid node
    const bool fluid;
    
    // whether a lattice node is a boundary node
    const bool boundary;
    
    // whether a lattice node is a solid node
    const bool solid;
    
    // macroscopic quantities
    // density of the fluid
    double rho;
    // velocity of the fluid
    std::array<double, dim> u;


public:
    LatticeNode(
        bool is_fluid_,
        bool is_boundary_,
        bool is_solid_
    ) :
    fluid (is_fluid_),
    boundary (is_boundary_),
    solid (is_solid_),
    populations ({0.0}),
    rho (0.0),
    rho_u ({0.0})
    {};
    
    ~LatticeNode() = default;

    const std::array<double, numVelocities> get_populations() const {
        return populations;
    }

    const std::array<double, dim> get_rho_u() const {
        return rho_u;
    }

    const double get_rho() const {
        return rho;
    }

    const bool is_fluid() const {
        return fluid;
    }

    const bool is_boundary() const {
        return boundary;
    }
    
    const bool is_solid() const {
        return solid;
    }

    /**
     * @param velocity_set the velocity set from which velocities are taken
     * @return the macroscopic density of the fluid and the velocity as a vector
    */
    const std::tuple<double, std::array<double, dim>> macroscopic_quantities(const VelocitySet& velocity_set) const {
        // calculating the density
        rho = 0.0;
        for (const auto& population : populations) 
        {
            rho += population;
        }

        u = {0.0};

        return [rho, u];
    }
    // TODO: setters
    // TODO: implement calculations for rho and rho_u
};

#endif  // HH_LBM_LATTICE_NODE