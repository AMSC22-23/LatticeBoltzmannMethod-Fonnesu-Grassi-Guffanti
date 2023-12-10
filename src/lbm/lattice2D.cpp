#include "lattice2D.hpp"

Lattice2D::Lattice2D(const std::string& input_dir_path, 
    const std::string& output_dir_path_, 
    const VelocitySet& velocity_set, 
    std::shared_ptr<CollisionModel> collision_model, 
    const double tau, 
    const double delta_t):
Lattice(input_dir_path, output_dir_path_, dim, velocity_set, collision_model, tau, delta_t)
{
    lattice_reader = std::make_unique<LatticeReader2D>(input_dir_path);
    if (!lattice_reader->read_lattice_structure(lattice, boundary_list , lattice_width, lattice_height))
    {
        std::cerr << "[!ERROR!] lattice structure matrix could not be read" << std::endl;
        assert(false);
    }
    initialize_lattice();
};

void Lattice2D::save_output_data(std::size_t iteration_count) const
{
    // file where the macroscopic density is saved
    std::string output_file_path_rho = output_dir_path + "/output-" + std::to_string(iteration_count) + "rho.txt" ;
    std::ofstream output_file_rho(output_file_path_rho);

    // file where the norm of the velocity is saved
    std::string output_file_path_u_x = output_dir_path + "/output-" + std::to_string(iteration_count) + "u_x.txt";
    std::ofstream output_file_u_x(output_file_path_u_x);

    // file where the norm of the velocity is saved
    std::string output_file_path_u_y = output_dir_path + "/output-" + std::to_string(iteration_count) + "u_y.txt";
    std::ofstream output_file_u_y(output_file_path_u_y);

    for (std::size_t i = 0; i < lattice_height; i++)
    {
        for (std::size_t j = 0; j < lattice_width; j++)
        {
            auto [rho, u] = lattice[i][j].get_macroscopic_quantities();
                        
            output_file_rho << rho << " ";
            output_file_u_x << u[0] << " ";
            output_file_u_y << u[1] << " ";

        }
        output_file_rho << std::endl;
        output_file_u_x << std::endl;
        output_file_u_y << std::endl;
        
    }
    
    output_file_rho.close();
    output_file_u_x.close();
    output_file_u_y.close();

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
    
    for (size_t i = 0; i < lattice_height; i++)
    {
        for (size_t j = 0; j < lattice_width; j++)
        {
            lattice[i][j].initialize_generic_boundary(velocity_set);
        }

    }
    
    if (!lattice_reader->read_lattice_input_velocities(lattice))
    {
        std::cout << "[!ERROR!] could not read input velocities" << std::endl;
        assert(false);
    }else if(!lattice_reader->read_lattice_input_rho(lattice)){
        std::cout << "[!ERROR!] could not read input densities" << std::endl;
        assert(false);
    }

    const WeightedDirection set_elements = velocity_set.get_velocity_set();

    const auto weights = set_elements.weight;

    // looping over all elements in the lattice
    for (std::size_t i = 0; i < lattice_height; ++i)
    {
        for (std::size_t j = 0; j < lattice_width; ++j)
        {
            if (lattice[i][j].is_fluid())
            {
                lattice[i][j].initialize_fluid_node(velocity_set, lattice[i][j].get_u(),lattice[i][j].get_rho());
            }
        }
    }

    std::cout << "LATTICE 2D:   lattice initialized" << std::endl;
}

void Lattice2D::perform_simulation_step() 
{
    // The simulation step is composed of many substeps:

    const double t_const = delta_t/tau;
    const double t_conj = 1-t_const;

    // 1. The equilibrium populations are calculated for each node
    // -> parallelizable!
    for (std::size_t i = 0; i < lattice_height; i++)
    {
        for (std::size_t j = 0; j < lattice_width; j++)
        {
            // TODO: can it be done only for fluid?
            if (lattice[i][j].is_fluid())
            {
                lattice[i][j].compute_equilibrium_populations(velocity_set.get_velocity_set());
            }

            // WHEN TO OUTPUT THE MACROSCOPIC QUANTITIES

            // 2. Perform the collisions
            if(lattice[i][j].is_fluid()){
                lattice[i][j].set_collision_populations() = collision_model->calc_collision(lattice[i][j].get_populations(), lattice[i][j].get_eq_populations(), t_const, t_conj);
            }
        }
    }

    // 3. Perform streaming
    perform_streaming();

    // 4. Perform the collision at the boundaries

    std::size_t size = boundary_list.size();
    for (std::size_t it = 0; it < size; it++)
    {  
        auto [i, j, type] = boundary_list[it];
        if (type == BOTTOM_WALL_2D) 
        {
            boundary_model.calc_bottom_wall_bounce(lattice[i][j]);
            lattice[i-1][j].set_population(2) = lattice[i][j].set_population(2);
            lattice[i-1][j-1].set_population(6) = lattice[i][j].set_population(6);
            lattice[i-1][j+1].set_population(5) = lattice[i][j].set_population(5);
        }
        else if (type == UPPER_WALL_2D )
        {
            boundary_model.calc_top_wall_bounce(lattice[i][j]);
            lattice[i+1][j].set_population(4) = lattice[i][j].set_population(4);
            lattice[i+1][j-1].set_population(7) = lattice[i][j].set_population(7);
            lattice[i+1][j+1].set_population(8) = lattice[i][j].set_population(8);
        }
        else if (type == LEFT_WALL_2D ) 
        {
            boundary_model.calc_left_wall_bounce(lattice[i][j]);
            lattice[i][j+1].set_population(1) = lattice[i][j].set_population(1);
            lattice[i-1][j+1].set_population(5) = lattice[i][j].set_population(5);
            lattice[i+1][j+1].set_population(8) = lattice[i][j].set_population(8);
        }
        else if (type == RIGHT_WALL_2D ) 
        {
            boundary_model.calc_right_wall_bounce(lattice[i][j]); 
            lattice[i][j-1].set_population(3) = lattice[i][j].set_population(3);
            lattice[i-1][j-1].set_population(6) = lattice[i][j].set_population(6);
            lattice[i+1][j-1].set_population(7) = lattice[i][j].set_population(7); 
        }
        else if (type == UPPER_LEFT_CORNER_2D ) 
        {
            boundary_model.calc_top_left_corner_bounce(lattice[i][j], lattice[i][j+1].get_rho());
            lattice[i][j+1].set_population(1) = lattice[i][j].set_population(1);
            lattice[i+1][j].set_population(4) = lattice[i][j].set_population(4);
            lattice[i+1][j+1].set_population(8) = lattice[i][j].set_population(8);
        }
        else if (type == UPPER_RIGHT_CORNER_2D ) 
        {
            boundary_model.calc_top_right_corner_bounce(lattice[i][j], lattice[i][j-1].get_rho()); 
            lattice[i][j-1].set_population(3) = lattice[i][j].set_population(3);
            lattice[i+1][j].set_population(4) = lattice[i][j].set_population(4);
            lattice[i+1][j-1].set_population(7) = lattice[i][j].set_population(7);
        }
        else if (type == BOTTOM_LEFT_CORNER_2D ) 
        {
            boundary_model.calc_bottom_left_corner_bounce(lattice[i][j], lattice[i][j+1].get_rho()); 
            lattice[i][j+1].set_population(1) = lattice[i][j].set_population(1);
            lattice[i-1][j].set_population(2) = lattice[i][j].set_population(2);
            lattice[i-1][j+1].set_population(5) = lattice[i][j].set_population(5);
        }
        else if (type == BOTTOM_RIGHT_CORNER_2D ) 
        {
            boundary_model.calc_bottom_right_corner_bounce(lattice[i][j], lattice[i][j-1].get_rho());
            lattice[i-1][j].set_population(2) = lattice[i][j].set_population(2);
            lattice[i][j-1].set_population(3) = lattice[i][j].set_population(3);
            lattice[i-1][j-1].set_population(6) = lattice[i][j].set_population(6);
        }
    }

    // 5. Update the macroscopic quantities
    for (std::size_t i = 0; i < lattice_height; i++)
    {
        for (std::size_t j = 0; j < lattice_width; j++)
        {
            if (lattice[i][j].is_fluid())
            {
                lattice[i][j].update_macroscopic_quantities(velocity_set);

            }
        }
    }
}

void Lattice2D::perform_streaming()
{
    for(std::size_t i = 0; i< lattice_height; i++)
    {
        for(std::size_t j = 0; j < lattice_width; j++)
        {
            // 3. Perform the streaming of fluid populations
            if (lattice[i][j].is_fluid()) {
                // streaming must be performed everywhere. In this way walls are able to perform bounces.
                lattice[i][j].set_population(0) = lattice[i][j].get_collision_populations()[0];
                lattice[i][j+1].set_population(1) = lattice[i][j].get_collision_populations()[1];
                lattice[i-1][j].set_population(2) = lattice[i][j].get_collision_populations()[2];
                lattice[i][j-1].set_population(3) = lattice[i][j].get_collision_populations()[3];
                lattice[i+1][j].set_population(4) = lattice[i][j].get_collision_populations()[4];
                lattice[i-1][j+1].set_population(5) = lattice[i][j].get_collision_populations()[5];
                lattice[i-1][j-1].set_population(6) = lattice[i][j].get_collision_populations()[6];
                lattice[i+1][j-1].set_population(7) = lattice[i][j].get_collision_populations()[7];
                lattice[i+1][j+1].set_population(8) = lattice[i][j].get_collision_populations()[8];
            }
        }
    }
}