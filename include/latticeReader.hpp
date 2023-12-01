#ifndef HH_LATTICE_READER
#define HH_LATTICE_READER

#include <iostream>
#include <string>
#include <filesystem>

/**
 * Abstract class that virtual methods that are overridden to implement the 
 * reading of a lattice structure and the values used to initialize the lattice
 * @author Luca Guffanti
*/
class LatticeReader
{
protected:
    /**
     * Path to the input directory which contains all the input files
    */
    const std::filesystem::path input_dir_path;
private:
    /**
     * Method that validates the path passed as input
    */
    bool validate_path();
public:
    LatticeReader(const std::string& input_dir_path_);
    virtual ~LatticeReader() = 0;
};



#endif // HH_LATTICE_READER