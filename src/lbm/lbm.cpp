#include "lbm.hpp"

lbm::lbm(std::size_t D, int Q, const std::string& input_dir_path, const std::string &collision_model, const std::string& boundary_model):
tau (1.0),
delta_t (1.0)
{
    velocity_set.initialize(D,Q);
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

    if(boundary_model == "BB")
    {
        boundary_ptr = std::make_shared<BounceBack>();
    }else if(collision_model == "NEBB")
    {
        boundary_ptr = std::make_shared<NEBB>();
    }else
    {
        std::cout << boundary_model << " not yet implemented" << std::endl;
    }

    switch (D){
        case 1:
            //TODO: to implement? punto di domanda?
            lattice_ptr=nullptr;
            break;
        case 2:
            lattice_ptr = std::make_unique<Lattice2D>(input_dir_path, output_dir_path, velocity_set, collision_ptr, tau, delta_t);
            break;
        case 3:
            //TODO: to implement
            lattice_ptr=nullptr;
            break;
        default:
            std::cout << "Errore" << std::endl;
            break;
    }
}

int lbm::compute(const std::size_t n_iter){

    for (std::size_t i = 0; i < n_iter; i++)
    {
        std::cout << "Iteration: " << i << " out of " << n_iter << ". ( " << i/(double)n_iter * 100.0 <<"% )" << std::endl;
        lattice_ptr->perform_simulation_step();
        lattice_ptr->save_output_data(i);
    }
    return 1;
    
}



