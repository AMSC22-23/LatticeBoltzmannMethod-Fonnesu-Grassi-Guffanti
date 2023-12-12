#include "lbm.hpp"

lbm::lbm(const std::size_t& D,
    const double reynolds_, 
    const std::string& input_dir_path, 
    const std::string& collision_model,
    const std::size_t& frequency_ = 1):
    frequency (frequency_),
    re (reynolds_)
    {
    std::string output_dir_path="../results";
    lattice_ptr=nullptr;

    if(collision_model == "BGK")
    {
        collision_ptr = std::make_shared<BGK>();
    }else if(collision_model == "TRT")
    {
        collision_ptr = std::make_shared<TRT>();
    }else if(collision_model == "MRT")
    {
        collision_ptr = std::make_shared<MRT>();
    }else
    {
        std::cout << collision_model << " not yet implemented" << std::endl;
    }

    switch (D)
    {
        case 2:
            velocity_set.initialize(2, 9);
            lattice_ptr = std::make_unique<Lattice2D>(input_dir_path, output_dir_path, velocity_set, collision_ptr, re);
            break;
        case 3:
            lattice_ptr=nullptr;
            break;
        default:
            std::cout << "[!ERROR!] only 2D and 3D are supported" << std::endl;
            break;
    }
}

void lbm::compute(const double time)
{
    double dt = re * (0.2 * 100 /re)/(100.0 * 100.0);
 
    const std::size_t n_iter = floor(time/dt);

    std::cout << "Time  : " << time << std::endl;
    std::cout << "dt    : " << dt << std::endl;
    std::cout << "iter  : " << n_iter << std::endl;

    std::size_t output_counter = 0;
    for (std::size_t i = 0; i <= n_iter; i++)
    {
        lattice_ptr->set_inlets(i);
        // REINITIALIZATION OF INLET FIELDS
        lattice_ptr->perform_simulation_step();
        if(i % frequency == 0)
        {
            std::cout << "iteration " << i << " out of " << n_iter << std::endl; 
            // print the output only if the iteration is a multiple of the frequency. 
            lattice_ptr->save_output_data(output_counter);
            output_counter++;
        }
    }    
}



