

#ifndef LLALBM_OPENBOUNDARYPOLICY_HPP
#define LLALBM_OPENBOUNDARYPOLICY_HPP

#include <vector>
#include <map>
#include <array>
#include <cassert>
#include <execution>
#include <openacc.h>
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
        Logger log;
    public:
        OpenBoundaryPolicy()
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
            /**
             * @brief This method does nothing as the open boundary condition does not require any update
             * 
             * @param populations population tensor (unused)
             * @param boundary_coord vector of coordinates of boundary nodes (unused)
             * @param global_rho density tensor (unused)
             * @param global_u velocity tensor (unused)
             */
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

            /**
             * @brief This method does nothing as the open boundary condition does not require any update
             * 
             * @param populations population tensor (unused)
             * @param boundary_coord vector of coordinates of boundary nodes (unused)
             * @param global_rho density tensor (unused)
             * @param global_u velocity tensor (unused)
             */
            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                return;
            }
    };
    
    template<std::size_t dim>
    class STDExecOpenBoundaryPolicy : public BoundaryPolicyTag, public STDExecTag
    {
    private:
        /* data */
        Logger log;
    public:
        STDExecOpenBoundaryPolicy(/* args */)
        : log("Generic Open Boundary" , std::cout)
        {
            log.error("Generic Open Boundary not implemented");
            assert(false);
        }
    };
    
    template<>
    class STDExecOpenBoundaryPolicy<2> : public BoundaryPolicyTag, public STDExecTag
    {       
        public: 

            /**
             * @brief This method does nothing as the open boundary condition does not require any update
             * 
             * @param populations population tensor (unused)
             * @param boundary_coord vector of coordinates of boundary nodes (unused)
             * @param global_rho density tensor (unused)
             * @param global_u velocity tensor (unused)
             */
            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                return;
            }
    };

    template<std::size_t dim>
    class OpenACCOpenBoundaryPolicy : public BoundaryPolicyTag, public OpenACCTag
    {
    private:
        Logger log;
    public:
        OpenACCOpenBoundaryPolicy()
        : log("Generic OpenBoundary" , std::cout)
        {
            log.error("Generic OpenBoundary not implemented");
            assert(false);
        }
    };
    
    template<>
    class OpenACCOpenBoundaryPolicy<2> : public BoundaryPolicyTag, public OpenACCTag
    {       
        public:
            /**
             * @brief This method does nothing as the open boundary condition does not require any update
             * 
             * @param populations population tensor (unused)
             * @param boundary_coord vector of coordinates of boundary nodes (unused)
             * @param global_rho density tensor (unused)
             * @param global_u velocity tensor (unused)
             */ 
            static void update_boundaries(Tensor<double, 3> &populations, std::vector<BoundaryPoint<2>> &boundary_coord, Tensor<double, 2> global_rho, Tensor<double, 3> global_u)
            {
                return;
            }
    };
    

} // namespace llalbm::core::boundaries

#endif // LLALBM_OpenBoundaryPolicy_HPP