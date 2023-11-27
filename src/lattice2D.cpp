#include "lattice2D.hpp"

Lattice2D::Lattice2D(const std::string& input_file_path_, const std::string& output_dir_path_, const VelocitySet& velocity_set):
Lattice(input_file_path_, output_dir_path_, dim, velocity_set)
{
    read_input_file();
    initialize_lattice();
};

void Lattice2D::save_output_data() const
{
    //TODO: not yet implemented
    std::cout << "LATTICE 2D:   saving output data" << std::endl;

}

void Lattice2D::read_input_file()
{
    std::cout << "LATTICE 2D:   reading input file " << input_file_path << std::endl;
    std::ifstream input_file(input_file_path);
    
    if (!input_file.is_open())
    {
        std::cerr << "LATTICE 2D:   could not open " << input_file_path << std::endl;
        assert(input_file.is_open());
    }

    std::string line;
    

    // Read the name of the lattice
    std::getline(input_file, lattice_name);
    std::cout << "  Lattice name: " << lattice_name << std::endl;

    // Read the dimensions of the lattice
    std::getline(input_file, line);
    std::size_t dim = std::stoul(line);
    if (dim != Lattice2D::dimensions)
    {
        std::cerr << "[ERROR]   wrong format (wrong dimensions)" << std::endl;
    }
    std::cout << "  Lattice dimensions: " << dim << std::endl;

    // Read the width of the lattice
    std::getline(input_file, line);
    std::size_t width = std::stoul(line);
    if (width == 0)
    {
        std::cerr << "[ERROR]   wrong format (no width)" << std::endl;
        assert(width > 0);
    }
    std::cout << "  Lattice width: " << width << std::endl;

    // Read the height of the lattice
    std::getline(input_file, line);
    std::size_t height = std::stoul(line);
    if (height == 0)
    {
        std::cerr << "[ERROR]   wrong format (no height)" << std::endl;
        assert(height > 0);
    }
    std::cout << "  Lattice height: " << height << std::endl;

    // Read data

    // first reserve the size of the vector
    lattice.resize(height, std::vector<LatticeNode<2>>(width));

    std::cout << "  reserved the size for the matrix" << std::endl;


    // the matrix is read from 
    Node_type type;
    int type_i;
    std::cout << "  reading data" << std::endl;

    for (std::size_t i = 0; i < height; i++)
    {
        for (std::size_t j = 0; j < width; j++)
        {
            // std::cout << "  element " << i << " " << j << std::endl;
            input_file >> type_i;
            type = static_cast<Node_type>(type_i);
            lattice[i][j].set_type() = type;
        }
    }

    std::cout << "  FINISHED READING!" << std::endl;
    input_file.close();

    lattice_height = height;
    lattice_width = width;
}

void Lattice2D::log_specific_data() const 
{
    //TODO: not yet implemented
    //std::cout << "  Lattice width  : " << Lattice2D::lattice.size() << std::endl;
    //std::cout << "  Lattice height : " << Lattice2D::lattice[0].size() << std::endl;
}

void Lattice2D::perform_simulation_step() 
{
    //TODO: not yet implemented
}

void Lattice2D::initialize_lattice() 
{
    std::cout << "LATTICE 2D:   initializing lattice" << std::endl;
    const std::vector<WeightedDirection> set_elements = velocity_set.get_velocity_set();

    for (std::size_t i = 0; i < lattice_height; ++i)
    {
        for (std::size_t j = 0; j < lattice_width; ++j)
        {
            // std::cout << i << " " << j << std::endl;
            lattice[i][j].initialize_node(set_elements);
        }
    }

    // looping over all elements in the lattice
    std::cout << "LATTICE 2D:   lattice initialized" << std::endl;
}