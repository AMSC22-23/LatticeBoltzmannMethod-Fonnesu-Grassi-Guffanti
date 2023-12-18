#include "latticeReader.hpp"

LatticeReader::LatticeReader(const std::string& input_dir_path_) :
input_dir_path (input_dir_path_),
input_lattice_path (input_dir_path_ + "/lattice.mtx"),
input_rho_path (input_dir_path_ + "/rho.mtx"),
input_u_path (input_dir_path_ + "/u")
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

    Eigen::MatrixXi mat(in);
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
                boundary_list[k] = {i, j, NONE};
                k++;
            }
        }
    }
    std::cout << "finished reading" << std::endl;

    // classifying the type of boundary
    std::cout << "classifying the type of boundary" << std::endl;

    // TODO: with the & operator the node is classified, but there are errors yet
    for(std::tuple<std::size_t, std::size_t, BoundaryType2DEnum> &boundary : boundary_list)
    {
        size_t i = std::get<0>(boundary);
        size_t j = std::get<1>(boundary);


        //CORNERS
        if(i!=0 && j!=0){
            if((lattice[i-1][j].is_generic_boundary() && lattice[i][j-1].is_generic_boundary()) &&
                !lattice[i-1][j-1].is_generic_boundary()){
                    std::get<2>(boundary) = BOTTOM_RIGHT_CORNER_2D;
                    // std::cout << i << " " << j << " " << std::get<2>(boundary) << std::endl;
                    continue;
                }
        }
        if(i!=(height-1) && j!=0){
            if((lattice[i+1][j].is_generic_boundary() && lattice[i][j-1].is_generic_boundary()) &&
                !lattice[i+1][j-1].is_generic_boundary()){
                    std::get<2>(boundary) = UPPER_RIGHT_CORNER_2D;
                    // std::cout << i << " " << j << " " << std::get<2>(boundary) << std::endl;
                    continue;
                }
        }
        if(i!=0 && j!=(width-1)){
            if((lattice[i-1][j].is_generic_boundary() && lattice[i][j+1].is_generic_boundary()) &&
                !lattice[i-1][j+1].is_generic_boundary()){
                    std::get<2>(boundary) = BOTTOM_LEFT_CORNER_2D;
                    // std::cout << i << " " << j << " " << std::get<2>(boundary) << std::endl;
                    continue;
                }
        }
        if(i!=(height-1) && j!=(width-1)){
            if((lattice[i+1][j].is_generic_boundary() && lattice[i][j+1].is_generic_boundary()) &&
                !lattice[i+1][j+1].is_generic_boundary()){
                    std::get<2>(boundary) = UPPER_LEFT_CORNER_2D;
                    // std::cout << i << " " << j << " " << std::get<2>(boundary) << std::endl;
                    continue;
                }
        }
        

        //WALLS
        if(i!=0 && i!=(height-1)){
            if(j!=0 ){
                if( (lattice[i+1][j].is_generic_boundary() && lattice[i-1][j].is_generic_boundary()) &&
                    (!lattice[i][j-1].is_generic_boundary())){
                        std::get<2>(boundary) = RIGHT_WALL_2D;
                        // std::cout << i << " " << j << " " << std::get<2>(boundary) << std::endl;
                        continue;
                }
            }
            if(j!=(width-1)){ // il meno uno c'è? punto di domanda FIXME:
                if( (lattice[i+1][j].is_generic_boundary() && lattice[i-1][j].is_generic_boundary()) &&
                    (!lattice[i][j+1].is_generic_boundary())){
                        std::get<2>(boundary) = LEFT_WALL_2D;
                        // std::cout << i << " " << j << " " << std::get<2>(boundary) << std::endl;
                        continue;
                }
            }
        }
        if(j!=0 && j!=(width-1)){
            if(i!=0){
                if( (lattice[i][j+1].is_generic_boundary() && lattice[i][j-1].is_generic_boundary()) &&
                    (!lattice[i-1][j].is_generic_boundary())){
                        std::get<2>(boundary) = BOTTOM_WALL_2D;
                        // std::cout << i << " " << j << " " << std::get<2>(boundary) << std::endl;
                        continue;
                }
            }
            if(i!=(height-1)){
                if( (lattice[i][j+1].is_generic_boundary() && lattice[i][j-1].is_generic_boundary()) &&
                    (!lattice[i+1][j].is_generic_boundary())){
                        std::get<2>(boundary) = UPPER_WALL_2D;
                        // std::cout << i << " " << j << " " << std::get<2>(boundary) << std::endl;
                        continue;
                }
            }
        }
    }
    std::cout << "boundary recognized" << std::endl;
    return true;
}

bool LatticeReader2D::read_lattice_input_rho(LatticeGrid2D& lattice)
{
    if (!(std::filesystem::exists(input_rho_path) && std::filesystem::is_regular_file(input_rho_path)))
    {
        std::cout << "density field was not found. Defaulting to 0.0 for boundaries and 1.0 for fluid" << std::endl;
        for (std::size_t i = 0; i < lattice.size(); ++i)
        {
            for (std::size_t j = 0; j < lattice[0].size(); ++j)
            {
                if (!lattice[i][j].is_generic_boundary())
                {
                    lattice[i][j].set_rho() = 1.0;
                } else 
                {
                    lattice[i][j].set_rho() = 0.0;
                }
            }
        }
    } else 
    {
        std::cout << "density field found" << std::endl;
        Eigen::SparseMatrix<double> in;
        if (!Eigen::loadMarket(in, input_rho_path))
        {
            std::cerr << "could not load density matrix file" << std::endl;
            assert(false);
        }

        Eigen::MatrixXd data(in);
        std::size_t rows = in.rows();
        std::size_t cols = in.cols();

        for (std::size_t i = 0; i < rows; ++i)
        {
            for (std::size_t j = 0; j < cols; ++j) 
            {
                lattice[i][j].set_rho() = data(i, j);
            }
        }
    }



    return true;
}

bool LatticeReader2D::read_lattice_input_velocities(LatticeGrid2D& lattice)
{
    const std::string u_x_path = input_u_path + "x.mtx";
    const std::string u_y_path = input_u_path + "y.mtx";


    bool x_present, y_present;
    // check existence of file
    if (!(x_present = std::filesystem::exists(u_x_path) && std::filesystem::is_regular_file(u_x_path)))
    {
        std::cout << "x-velocity field was not found. Defaulted to 0 everywhere" << std::endl;
    } else
    {
        std::cout << "x-velocity field found" << std::endl;
    }
    if (!(y_present = std::filesystem::exists(u_y_path) && std::filesystem::is_regular_file(u_y_path)))
    {
        std::cout << "y-velocity field was not found. Defaulted to 0 everywhere" << std::endl;
    } else
    {
        std::cout << "y-velocity field found" << std::endl;
    }

    if (!x_present)
    {
        
        std::cout << "Defaulting x velocity to 0" << std::endl;
        for (std::size_t i = 0; i < lattice.size(); ++i)
        {
            for (std::size_t j = 0; j < lattice[0].size(); ++j)
            {
                lattice[i][j].set_u(0) = 0.0;
            }
        }
    }
    if (!y_present)
    {
        
        std::cout << "Defaulting y velocity to 0" << std::endl;
        for (std::size_t i = 0; i < lattice.size(); ++i)
        {
            for (std::size_t j = 0; j < lattice[0].size(); ++j)
            {
                lattice[i][j].set_u(1) = 0.0;
            }
        }
    }

    if (x_present)
    {
        if (!load_2D_velocity_matrix(lattice, u_x_path, 0))
        {
            std::cerr << "x-velocity input matrix could not be read" << std::endl;
            assert(false);
        }
    } 
    if (y_present)
    {
        if (!load_2D_velocity_matrix(lattice, u_y_path, 1))
        {
            std::cerr << "y-velocity input matrix could not be read" << std::endl;
            assert(false);
        }
    }
    return true;
}

bool LatticeReader2D::load_2D_velocity_matrix(LatticeGrid2D& lattice, const std::string& path, const std::size_t& dir)
{
    Eigen::SparseMatrix<double> mat;
    if (!Eigen::loadMarket(mat, path))
    {
        return false;
    }
    Eigen::MatrixXd mat_dense(mat);

    std::size_t rows = mat.rows();
    std::size_t cols = mat.cols();
    for (std::size_t i = 0; i < rows; ++i)
    {
        for (std::size_t j = 0; j < cols; ++j)
        {
            lattice[i][j].set_u(dir) = mat_dense(i, j);
        }
    }
    return true;
}
