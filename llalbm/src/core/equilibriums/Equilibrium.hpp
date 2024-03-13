// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
// ======================================

// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "../../utils/loggers/Logger.hpp"
#include "../../utils/aliases.hpp"
// =======================================

namespace llalbm::core::equilibrium
{

    template<
        std::size_t dim>
    class Equilibrium
    {};

    template<>
    class Equilibrium<2>{
        private:
        double rho,ux,uy,u2;
        
        public:

        void calc_equilibrium(Tensor<double, 3> &populations, Tensor<double, 3> &equilibrium_populations, Tensor<double, 3> &global_u, Tensor<double,2> &global_rho){
            constexpr double two9 = (2.0/9.0);
            constexpr double one18 = (1.0/18.0);
            constexpr double one36 = (1.0/36.0);

            for (size_t i = 0; i < populations.dimension(0); i++)
            {
                for (size_t j = 0; j < populations.dimension(1); j++)
                {   
                    rho = global_rho(i,j);
                    ux = global_u(i,j,0);
                    uy = global_u(i,j,1);
                    u2 = ux * ux + uy * uy;
                    equilibrium_populations(i,j,0) =  two9 * rho * (2 - 3*u2);
                    equilibrium_populations(i,j,1) =  one18 * rho * (2 + 6 * ux + 9 * ux * ux - 3 * u2);
                    equilibrium_populations(i,j,2) =  one18 * rho * (2 + 6 * uy + 9 * uy * uy - 3 * u2);
                    equilibrium_populations(i,j,3) =  one18 * rho * (2 - 6 * ux + 9 * ux * ux - 3 * u2);
                    equilibrium_populations(i,j,4) =  one18 * rho * (2 - 6 * uy + 9 * uy * uy - 3 * u2);
                    equilibrium_populations(i,j,5) = one36 * rho * (1 + 3 * (ux + uy) + 9*ux*uy + 3*u2);
                    equilibrium_populations(i,j,6) = one36 * rho * (1 - 3 * (ux - uy) - 9*ux*uy + 3*u2);
                    equilibrium_populations(i,j,7) = one36 * rho * (1 - 3 * (ux + uy) + 9*ux*uy + 3*u2);
                    equilibrium_populations(i,j,8) = one36 * rho * (1 + 3 * (ux - uy) - 9*ux*uy + 3*u2);
                }
                
            }
            
        }
    };
} // namespace llalbm::core::equilibrium