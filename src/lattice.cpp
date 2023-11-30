#include "lattice.hpp"

Lattice::Lattice(const std::string& input_file_path_,
 const std::string& output_dir_path_, const int dimensions_, const VelocitySet& velocity_set_) :
input_file_path (input_file_path_),
output_dir_path (output_dir_path_),
dimensions (dimensions_),
velocity_set (velocity_set_)
{
    create_output_directory();
};

void Lattice::log_data() const
{
    std::cout << "LATTICE " << std::endl;
    std::cout << "  Dimensions  : " << dimensions << std::endl;
    std::cout << "  Input file  : " << input_file_path << std::endl;
    std::cout << "  Output dir  : " << output_dir_path << std::endl;
    log_specific_data();
}

void Lattice::create_output_directory() 
{
    std::filesystem::path path;
    // if the path is not correct than the program defaults to ./results
    if (std::filesystem::status(output_dir_path).type() != std::filesystem::file_type::directory)
    {
        path = "../results/";
        output_dir_path = "../results";
        std::cout << "[WARNING] path to output directory was invalid. Defaulting to ../results" << std::endl;
        // and if the directory already exists it isn't created
        if (!std::filesystem::exists(path))
        {
            bool status;
            if (!(status = std::filesystem::create_directory(path)))
            {
                std::cerr << "[!ERROR!]   could not create results directory" << std::endl;
                assert(status == true);
            }
        } else
        {
            std::cout << "-> results directory already exists" << std::endl;
        }
    }
}