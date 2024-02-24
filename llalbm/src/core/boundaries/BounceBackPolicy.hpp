/**
 * @file BounceBackPolicy.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-02-20
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <vector>
#include <map>
#include <array>
// =========== EIGEN INCLUDES ===========
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <Eigen/Core>
// ======================================

// =========== LLALBM INCLUDES ===========
#include "Logger.hpp"
// =======================================


namespace llalbm::core::boundaries
{
    
    using namespace llalbm::util::logger;
    template<std::size_t d>
    using Point = Matrix<std::size_t, d, 1>;

    template<
        std::size_t dim>
    class BounceBackPolicy
    {
    private:
        /* data */
    public:
        BounceBackPolicy(/* args */)
        : log("Generic BB" , std::cout)
        {
            log.error("Generic BB not implemented");
            assert(false);
        }
    };

    template<>
    class BounceBackPolicy<2>{
        
        private:
            std::array<Eigen::Index, 2> lattice_nodes; 
            std::vector<Point<2>> boundary_nodes;
            
        public: 
            void load_nodes(std::array<Eigen::Index, 2> &l, std::vector<Point<2>> &b){
                lattice_nodes = l;
                boundary_nodes = b;
            }

            void update_boundaries(){
                std::cout << "update boundaries to do" << std::endl;
            }
    };



} // namespace llalbm::core::boundaries
