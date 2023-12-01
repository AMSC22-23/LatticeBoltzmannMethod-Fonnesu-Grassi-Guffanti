#ifndef HH_LATTICE_READER_2D
#define HH_LATTICE_READER_2D

#include "latticeReader.hpp"
#include <filesystem>
#include <iostream>

/**
 * Implementation of the abstract LatticeReader class
 * used to read the sets of matrices that describe a 2D lattice and
 * the initialization values
*/
class LatticeReader2D : public LatticeReader
{
private:

public:
    /**
     * Constructs the object
     * @param input_dir_path the path to the input directory
    */
    LatticeReader2D(const std::string input_dir_path_);
    virtual ~LatticeReader2D() = default;
};



#endif // HH_LATTICE_READER_2D