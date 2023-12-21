#include <iostream>
#include "lbm.hpp"

int main(int argc , char** argv)
{
    const int required_arguments = 4;
    if(argc < required_arguments){
        std::cout<< "Arguments :  dimensions path_to_input_to_dir collision_model [-r reynolds] [-f frequency to print data]" << std::endl;
        return 0;
    }
    double r=100.0;
    int frequency = 100;

    int d= std::stoi(argv[1]);
    std::string input_dir = argv[2];
    std::string collision_model = argv[3];

    if(argc > required_arguments){
        for(int i= required_arguments; i < argc ; i++){
            if(std::strcmp(argv[i],"-r") == 0){
                try{
                    r = std::stod(argv[++i]);
                }
                catch(const std::logic_error& ignored){};
            }else if(std::strcmp(argv[i],"-f") == 0){
                try{
                    frequency = std::stoi(argv[++i]);
                }
                catch(const std::logic_error& ignored){};
            }
        }
    }
    lbm LBM(d, r, input_dir,collision_model, frequency);
    LBM.compute(20.0);
    return 0;
}
