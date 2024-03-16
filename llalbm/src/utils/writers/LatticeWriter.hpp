#ifndef LLALBM_LATTICEWRITER_HPP
#define LLALBM_LATTICEWRITER_HPP

// =========== STL INCLUDES ===========
#include <filesystem>
#include <cassert>
#include <fstream>
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

namespace llalbm::util::writer
{
    using namespace llalbm::util::logger;

    static Logger logger("LatticeWriter", std::cout);
    static std::string output_dir_path = "../results/";

    
    void create_output_directory()
    {
        std::filesystem::path path;

        path = output_dir_path;

        if (!std::filesystem::exists(path))
        {
            assert(std::filesystem::create_directory(path) && "ERROR: could not create results directory");
            logger.info("  -> results directory built at path " + output_dir_path);
        }else
        {
            logger.info("-> results directory already exists");
            logger.info("clearing the directory");

            for (const auto& entry : std::filesystem::directory_iterator(output_dir_path))
            {
                std::cout << entry.path();
                if (entry.path().extension() == ".bin" || entry.path().extension() == ".txt")
                {
                    std::filesystem::remove(entry.path());
                    std::cout<< " removed" << std::endl;
                }
            }
        }
    }

    template<std::size_t dim>
    void write_lattice_file(
        const Tensor<double, dim + 1>& global_u, 
        std::size_t iteration_count)
    {
        std::string output_file_path_u_x = output_dir_path + "/output-" + std::to_string(iteration_count) + "u_x.txt";
        std::ofstream output_file_u_x(output_file_path_u_x);

        std::string output_file_path_u_y = output_dir_path + "/output-" + std::to_string(iteration_count) + "u_y.txt";
        std::ofstream output_file_u_y(output_file_path_u_y);
        

        assert(output_file_u_x.is_open() && "ERROR: u_x file could not be opened");
        assert(output_file_u_y.is_open() && "ERROR: u_y file could not be opened");

        for (Eigen::Index i = 0; i < global_u.dimension(0); ++i) {
            for (Eigen::Index j = 0; j < global_u.dimension(1); ++j) {
                output_file_u_x << global_u(i, j, 0) << " ";
                output_file_u_y << global_u(i, j, 1) << " ";

            }
            output_file_u_x << std::endl;
            output_file_u_y << std::endl;

        }

        output_file_u_x.close();
        output_file_u_y.close();

        logger.info("Iteration " + std::to_string(iteration_count));
    }

}

#endif