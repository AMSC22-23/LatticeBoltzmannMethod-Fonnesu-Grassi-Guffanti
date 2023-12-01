#include "lattice2D.hpp"

Lattice2D::Lattice2D(const std::string& input_dir_path, const std::string& output_dir_path_, const VelocitySet& velocity_set, std::shared_ptr<CollisionModel> collision_model, std::shared_ptr<Boundary> boundary_model, const double tau, const double delta_t):
lattice_reader (input_dir_path),
Lattice(input_dir_path, output_dir_path_, dim, velocity_set, collision_model, boundary_model, tau, delta_t)
{
    lattice_reader.read_lattice_structure();
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
            auto [rho, u] = lattice(i, j).get_macroscopic_quantities();
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
            if (lattice(i, j).is_fluid())
            {
                lattice(i, j).initialize_fluid_node(weights, {0.0, 0.0}, 1.0);
            } else if (lattice(i, j).is_open_boundary())
            {
            } else 
            {
                lattice(i, j).initialize_generic_boundary(velocity_set);
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
            lattice(i, j).compute_equilibrium_populations(velocity_set.get_velocity_set());

            // WHEN TO OUTPUT THE MACROSCOPIC QUANTITIES

            // 2. Perform the collisions
            if(lattice(i, j).is_fluid()){
                // lattice(i, j).set_collision_populations() = collision_model->calc_collision();
            }

            // 3. Perform the streaming
            perform_streaming();

            // 4. Perform the propagation at the boundaries

            // 5. Update the macroscopic quantities
            lattice(i, j).update_macroscopic_quantities(velocity_set);

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