#include "lbm.hpp"

lbm::lbm(std::size_t D, int Q, const std::string& input_file_path, const std::string &collision_model, const std::string& boundary_model):
tau (1.0),
delta_t (1.0)
{
    velocity_set.initialize(D,Q);
    std::string output_dir_path="../output";
    lattice_ptr=nullptr;
    switch (D){
        case 1:
            //TODO: to implement? punto di domanda?
            lattice_ptr=nullptr;
            break;
        case 2:
            lattice_ptr = std::make_unique<Lattice2D>(input_file_path, output_dir_path,velocity_set);
            break;
        case 3:
            //TODO: to implement
            lattice_ptr=nullptr;
            break;
        default:
            std::cout << "Errore" << std::endl;
            break;
    }

    /*switch (collision_model){
        case "BGK":
            collision_ptr = std::make_unique<BGK>();
            break;
        case "TRT":
            collision_ptr = std::make_unique<TRT>();
            break;
        case "MRT":
            collision_ptr = std::make_unique<MRT>();
            break;
        default: std::cout << collision_model << " not yet implemented" << std::endl;
    }
    switch (boundary_model){
    case "BB":
        boundary_ptr = std::make_unique<BounceBack>();
        break;
    
    default:
        std::cout << boundary_model << " not yet implemented" << std::endl;
        break;
    }*/

    //lattice_ptr->initialize_lattice();
}

int lbm::compute(const int n_iter){

    for (size_t i = 0; i < n_iter; i++)
    {
        lattice_ptr->perform_simulation_step();
        //lattice_ptr->save_output_data();
    }
    return 1;
    
}

