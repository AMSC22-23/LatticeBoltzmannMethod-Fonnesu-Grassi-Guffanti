#include "lbm.hpp"

lbm::lbm(size_t D, int Q, const std::string& input_file_path):
tao (1.0),
deltaT (1.0)
{
    velocitySet.initialize(D,Q);
    std::string output_dir_path="../output";
    lattice_ptr=nullptr;
    switch (D){
        case 1:
            //TODO: to implement? punto di domanda?
            lattice_ptr=nullptr;
            break;
        case 2:
            lattice_ptr = std::make_unique<Lattice2D>(input_file_path, output_dir_path);
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

