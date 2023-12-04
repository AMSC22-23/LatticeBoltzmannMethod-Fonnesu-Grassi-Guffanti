#include "latticeReader.hpp"

LatticeReader::LatticeReader(const std::string& input_dir_path_) :
input_dir_path (input_dir_path_),
input_lattice_path (input_dir_path_ + "/lattice.mtx"),
input_rho_path (input_dir_path_ + "/rho.mtx"),
input_u_path (input_dir_path_ + "/u.mtx")
{   
    if (!validate_path())
    {
        std::cout << "[!ERROR!] Path to input directory is invalid" << std::endl;
        assert(false);
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

bool LatticeReader2D::read_lattice_structure(LatticeGrid2D& lattice, 
    BoundaryList2D& boundary_list,
    std::size_t& width, 
    std::size_t& height)
{
    Eigen::SparseMatrix<int> in;

    // loading the matrix in the mtx format
    std::cout << "reading data" << std::endl;
    if (!Eigen::loadMarket(in, input_lattice_path))
    {  
        return false;
    }

    width = in.cols();
    height = in.rows();
    std::cout << "width of the domain      : " << width << std::endl;
    std::cout << "height of the domain     : " << height << std::endl;
    std::cout << "number of boundary nodes : " << in.nonZeros() << "(" << in.nonZeros()/(double)in.size() * 100.0 << "%)" << std::endl;

    lattice.resize(in.rows(), std::vector<LatticeNode<2>>(in.cols()));
    boundary_list.resize(in.nonZeros());

    Node_type type; 
    int type_i;

    Eigen::MatrixX<int> mat(in);
    std::size_t k = 0;
    for (std::size_t i = 0; i < height; ++i)
    {
        for (std::size_t j = 0; j < width; ++j)
        { 
            type_i = in.coeffRef(i, j);
            type = static_cast<Node_type>(type_i);
            lattice[i][j].set_type() = type;
            if (type != FLUID)
            {
                boundary_list[k] = {lattice[i][j], i, j};
                k++;
            }
        }
    }
    std::cout << "finished reading" << std::endl;
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