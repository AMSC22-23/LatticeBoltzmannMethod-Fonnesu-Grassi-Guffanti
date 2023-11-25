#include "lattice.hpp"

Lattice::Lattice(const std::string& input_file_path_, const std::string& output_dir_path_, const int dimensions_) :
input_file_path (input_file_path_),
output_dir_path (output_dir_path_),
dimensions (dimensions_)
{};

void Lattice::log_data() const
{
    std::cout << "LATTICE " << std::endl;
    std::cout << "  Dimensions  : " << dimensions << std::endl;
    std::cout << "  Input file  : " << input_file_path << std::endl;
    std::cout << "  Output dir  : " << output_dir_path << std::endl;
    log_specific_data();
}