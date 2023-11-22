#include <string>
#include <fstream>
#include "velocitySet.hpp"


VelocitySet::VelocitySet(const std::size_t D_, const int Q_):
    D(D_),
    Q(Q_)
    {
        std::string type="\"D" + std::to_string(D) + "Q" + std::to_string(Q)+"\": {";
        std::ifstream Sets("../resources/velocitySets.json");


        std::cout<< type;
        if(!Sets.is_open()){
            std::cout<< "Velocity sets not found." << std::endl;
        }else{
            std::string line;
            while (std::getline(Sets, line)) {
                if(line == type){
                    std::cout<< line << std::endl;
                }
                if (Sets.eof()) {
                    std::cout << "Hai letto tutto il file." << std::endl;
                    break;
                }
            }
            
            
        }
    };
    


