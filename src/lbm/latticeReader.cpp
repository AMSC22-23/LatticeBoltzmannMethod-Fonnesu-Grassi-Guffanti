#include "latticeReader.hpp"

LatticeReader::LatticeReader(const std::string& input_dir_path_) :
input_dir_path (input_dir_path_),
input_lattice_path (input_dir_path_ + "/lattice.mtx"),
input_rho_path (input_dir_path_ + "/rho.mtx"),
input_u_path (input_dir_path_ + "/u.mtx")
{   
    bool result;
    if (!(result = validate_path()))
    {
        std::cout << "[!ERROR!] Path to input directory is invalid" << std::endl;
        assert(result == true);
    }

    std::cout << "LatticeReader2D   correctly initialized" << std::endl;
    std::cout << "      input lattice structure path : " << input_lattice_path << std::endl;
    std::cout << "      input density field path     : " << input_rho_path << std::endl;
    std::cout << "      input velocity fields path   : " << input_u_path << std::endl;
}

bool LatticeReader::validate_path()
{
    return std::filesystem::exists(input_dir_path) &&
    std::filesystem::is_directory(input_dir_path);    
}

LatticeReader2D::LatticeReader2D(const std::string& input_dir_path_):
LatticeReader(input_dir_path_)
{};

bool LatticeReader2D::read_lattice_structure(Eigen::Matrix<LatticeNode<2>, Eigen::Dynamic, Eigen::Dynamic>& lattice, std::size_t& width, std::size_t& height)
{
    Eigen::SparseMatrix<int> in;

    // loading the matrix in the mtx format
    std::cout << "reading data" << std::endl;
    if (Eigen::loadMarket(in, input_lattice_path))
    {  
        return false;
    }

    width = in.cols();
    height = in.rows();
    

    Node_type type; 
    int type_i;
    Eigen::Matrix<int, -1, -1> dense_in = in; 
    for (std::size_t i = 0; i < height; ++i)
    {
        for (std::size_t j = 0; j < width; ++i)
        {
            type_i = dense_in(i, j);
            type = static_cast<Node_type>(type_i);
            lattice(i, j).set_type() = type;
        }
    }
    return true;
}

bool LatticeReader2D::read_lattice_input_rho()
{
    return true;
}

bool LatticeReader2D::read_lattice_input_velocities()
{
    return true;
}