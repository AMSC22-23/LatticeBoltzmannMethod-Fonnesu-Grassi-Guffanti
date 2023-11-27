#include <string>
#include <fstream>
#include "velocitySet.hpp"

void VelocitySet::set_velocity_set()
    {
        std::string type="D" + std::to_string(D) + "Q" + std::to_string(Q);
        std::ifstream Sets("../resources/velocitySets");

        if(!Sets.is_open()){
            std::cout<< "Velocity sets not found." << std::endl;
        }else{
            std::string line;
            while (std::getline(Sets, line)) {
                if(line == type){

                    //Directions reading

                    std::getline(Sets, line);
                    size_t separator =line.find(';');
                    while(line.length()!=0){
                        std::string direction=line.substr(0,separator);
                        line=line.substr(separator+1,line.length());

                        std::vector<int> dir;
                        while(direction.length()!=0){
                            size_t temp=direction.find(',');
                            
                            dir.push_back(std::stoi(direction.substr(0,temp)));
                            direction=direction.substr(temp+1,direction.length());
                            //std::cout << direction << " " << std::endl;
                        }
                        WeightedDirection newDir;
                        newDir.direction=dir;
                        newDir.weight=1;
                        this->Set.push_back(newDir);
                        separator =line.find(';');
                    }

                    //Weight reading
                    
                    std::getline(Sets, line);
                    separator=line.find(",");
                    for (auto &&currDir : this->Set)
                    {
                        
                        std::string weight=line.substr(0,separator);
                        line=line.substr(separator+1,line.length());
                        
                        //std::cout << weight << std::endl;

                        size_t temp=weight.find('/');
                        double numerator= std::stod(weight.substr(0,temp));
                        double denominator= std::stod(weight.substr(temp+1,weight.length()));
                        
                        double res=numerator/denominator;
                        //std::cout<< numerator << "/" << denominator << "=" << numerator/denominator << std::endl;
                        currDir.weight = res;
                        
                        separator=line.find(",");
                        
                    }
                    
                    
                    /* STAMPA VETTORE DIREZIONE E PESI   
                    for (auto &&i : this->Set){   
                        for (auto &&j : i.direction){
                            std::cout<< j << " ";
                        }
                        std::cout <<" "<< i.weight <<std::endl<< std::endl;
                    }*/
                }
                if (Sets.eof()) {
                    std::cout << "Hai letto tutto il file." << std::endl;
                    break;
                }
            }
            
            
        }
    };
    


