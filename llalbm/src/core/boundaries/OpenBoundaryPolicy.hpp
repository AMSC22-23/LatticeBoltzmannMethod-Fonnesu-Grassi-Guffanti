/**
 * @file OpenBoundaryPolicy.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-02-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef LLALBM_OPENBOUNDARYPOLICY_HPP
#define LLALBM_OPENBOUNDARYPOLICY_HPP

#include <vector>
#include <map>
#include <array>
#include <cassert>
// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "../../utils/loggers/Logger.hpp"
#include "../../utils/aliases.hpp" 
#include "../PolicyTypes.hpp"
// =======================================


namespace llalbm::core::boundaries
{
    using namespace llalbm::util::logger;
    template<std::size_t d>
    using Point = Matrix<std::size_t, d, 1>;

    template<std::size_t dim>
    class OpenBoundaryPolicy : public BoundaryPolicyTag, public SequentialTag
    {
    private:
        /* data */
        Logger log;
    public:
        OpenBoundaryPolicy(/* args */)
        : log("Generic OpenBoundary" , std::cout)
        {
            log.error("Generic OpenBoundary not implemented");
            assert(false);
        }
    };
    
    template<>
    class OpenBoundaryPolicy<2> : public BoundaryPolicyTag, public SequentialTag
    {       
        public: 
            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                return;
            }
    };
    

    template<std::size_t dim>
    class OMPOpenBoundaryPolicy : public BoundaryPolicyTag, public OMPTag
    {
    private:
        /* data */
        Logger log;
    public:
        OMPOpenBoundaryPolicy(/* args */)
        : log("Generic Open Boundary" , std::cout)
        {
            log.error("Generic Open Boundary not implemented");
            assert(false);
        }
    };
    
    template<>
    class OMPOpenBoundaryPolicy<2> : public BoundaryPolicyTag, public OMPTag
    {       
        public: 

            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                return;
            }
    };
    
} // namespace llalbm::core::boundaries

#endif // LLALBM_OpenBoundaryPolicy_HPP