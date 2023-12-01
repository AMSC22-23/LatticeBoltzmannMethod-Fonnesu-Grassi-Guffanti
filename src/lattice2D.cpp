#include "lattice2D.hpp"

Lattice2D::Lattice2D(const std::string& input_file_path_, const std::string& output_dir_path_, const VelocitySet& velocity_set, std::shared_ptr<CollisionModel> collision_model, std::shared_ptr<Boundary> boundary_model, const double tau, const double delta_t):
Lattice(input_file_path_, output_dir_path_, dim, velocity_set, collision_model, boundary_model, tau, delta_t)
{
    read_input_file();
    initialize_lattice();
};

void Lattice2D::save_output_data(std::size_t iteration_count) const
{
    // file where the macroscopic density is saved
    std::string output_file_path_rho = output_dir_path + "/output-" + std::to_string(iteration_count) + "rho.txt" ;
    std::ofstream output_file_rho(output_file_path_rho);

    // file where the norm of the velocity is saved
    std::string output_file_path_u = output_dir_path + "/output-" + std::to_string(iteration_count) + "u.txt";
    std::ofstream output_file_u(output_file_path_u);

    for (std::size_t i = 0; i < lattice_height; i++)
    {
        for (std::size_t j = 0; j < lattice_width; j++)
        {
            auto [rho, u] = lattice[i][j].get_macroscopic_quantities();
            // computing the two 2-norm of the vector
            auto two_norm = std::sqrt(
                u[0]*u[0] + u[1]*u[1]
            );

            output_file_rho << rho << " ";
            output_file_u << two_norm << " ";
        }
        output_file_rho << std::endl;
        output_file_u << std::endl;
    }

    output_file_rho.close();
    output_file_u.close();
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

    // Read the width of the lattice
    input_file >> line;
    std::size_t width = std::stoul(line);
    if (width == 0)
    {
        std::cerr << "[ERROR]   wrong format (no width)" << std::endl;
        assert(width > 0);
    }
    std::cout << "  Lattice width: " << width << std::endl;
    
    // Read the height of the lattice
    input_file >> line;
    std::size_t height = std::stoul(line);
    if (height == 0)
    {
        std::cerr << "[ERROR]   wrong format (no height)" << std::endl;
        assert(height > 0);
    }
    std::cout << "  Lattice height: " << height << std::endl;

    // Read the number of non fluid nodes in the domain
    input_file >> line;
    std::size_t nnz = std::stoul(line);
    std::cout << "  Non fluid nodes: " << line << "(" << (double)nnz/((double)width*height) * 100.0 << "%)" <<std::endl;
    // Read data

    // first reserve the size of the vector
    lattice.resize(height, std::vector<LatticeNode<2>>(width));
    lattice_height = height;
    lattice_width = width;

    std::cout << "  Reserved the size for the matrix" << std::endl;

    // the matrix is read from 
    Node_type type;
    std::cout << "  Reading data" << std::endl;

    // for each non zero element in the input matrix
    std::size_t x, y;
    std::size_t type_i;
    for (std::size_t i = 0; i < nnz; i++)
    {
        // read the x coordinate
        input_file >> line;
        x = std::stoul(line);
        // read the y coordinate
        input_file >> line;
        y = std::stoul(line);
        // read the type
        input_file >> line;
        type_i = std::stoul(line);
        
        type = static_cast<Node_type>(type_i);
        lattice[y][x].set_type() = type;
    }

    std::cout << "  FINISHED READING!" << std::endl;
    input_file.close();


}

void Lattice2D::log_specific_data() const 
{
    std::cout << "  Lattice Name   : " << lattice_name << std::endl;
    std::cout << "  Lattice width  : " << lattice_width << std::endl;
    std::cout << "  Lattice height : " << lattice_height << std::endl;
}

void Lattice2D::initialize_lattice() 
{
    std::cout << "LATTICE 2D:   initializing lattice" << std::endl;
    std::cout << "              reading input data from file" << std::endl;

    const WeightedDirection set_elements = velocity_set.get_velocity_set();

    const auto weights = set_elements.weight;

    // looping over all elements in the lattice
    for (std::size_t i = 0; i < lattice_height; ++i)
    {
        for (std::size_t j = 0; j < lattice_width; ++j)
        {
            if (lattice[i][j].is_fluid())
            {
                lattice[i][j].initialize_fluid_node(weights);
            } else if (lattice[i][j].is_open_boundary())
            {
            } else 
            {
                lattice[i][j].initialize_generic_boundary();
            }
        }
    }

    std::cout << "LATTICE 2D:   lattice initialized" << std::endl;
}

void Lattice2D::perform_simulation_step() 
{
    // The simulation step is composed of many substeps:

    // 1. The equilibrium populations are calculated for each node
    // -> parallelizable!
    for (std::size_t i = 0; i < lattice_height; i++)
    {
        for (std::size_t j = 0; j < lattice_width; j++)
        {
            lattice[i][j].compute_equilibrium_populations(velocity_set.get_velocity_set());

            // WHEN TO OUTPUT THE MACROSCOPIC QUANTITIES

            // 2. Perform the collisions
            if(lattice[i][j].is_fluid()){
                perform_collisions();
                //std::vector<double> fTemp = collision_model->calc_collision();
            }

            // 3. Perform the streaming
            perform_streaming();

            // 4. Perform the propagation at the boundaries

            // 5. Update the macroscopic quantities
            lattice[i][j].update_macroscopic_quantities(velocity_set);

        }
    }
}

void Lattice2D::perform_collisions() 
{
    // TODO: not yet implemented
}

void Lattice2D::perform_streaming()
{
    // TODO: not yet implemented
}