#include "lattice.hpp"

Lattice::Lattice(const std::string& input_dir_path_,
const std::string& output_dir_path_, const int dimensions_, const VelocitySet& velocity_set_, std::shared_ptr<CollisionModel> collision_model_, const double reynolds_) :
input_dir_path (input_dir_path_),
output_dir_path (output_dir_path_),
dimensions (dimensions_),
velocity_set (velocity_set_), 
collision_model (collision_model_),
reynolds (reynolds_)
{
    create_output_directory();
};

void Lattice::log_data() const
{
    std::cout << "LATTICE " << std::endl;
    std::cout << "  Dimensions  : " << dimensions << std::endl;
    std::cout << "  Input dir   : " << input_dir_path << std::endl;
    std::cout << "  Output dir  : " << output_dir_path << std::endl;
    log_specific_data();
}

void Lattice::create_output_directory() 
{
    std::filesystem::path path;
    std::cout << "trying to move to inserted input directory to produce output" << std::endl;
    std::cout << "  " << input_dir_path << std::endl;

    // if the path is not correct than the program defaults to ./results
    if (std::filesystem::is_directory(input_dir_path))
    {
        output_dir_path = input_dir_path + "/results/";
    } else {
        output_dir_path = "../results";
        std::cout << "[WARNING] path to output directory was invalid. Defaulting to ../results" << std::endl;
    }
    path = output_dir_path;
    // and if the directory already exists it isn't created

    if (!std::filesystem::exists(path))
    {
        bool status;
        if (!(status = std::filesystem::create_directory(path)))
        {
            std::cerr << "[!ERROR!]   could not create results directory" << std::endl;
            assert(status == true);
        }
        std::cout << "  -> results directory built at path " << output_dir_path << std::endl;
    } else
    {
        std::cout << "-> results directory already exists" << std::endl;
        std::cout << "clearing the directory (only .txt files)" << std::endl;

        for (const auto& entry : std::filesystem::directory_iterator(output_dir_path))
        {
            std::cout << entry.path();
            if (entry.path().extension() == ".txt")
            {
                std::filesystem::remove(entry.path());
                std::cout<< " removed" << std::endl;
            }
        }

    }
}

std::size_t& Lattice::set_omp_num_threads()
{
    return omp_num_threads;
}