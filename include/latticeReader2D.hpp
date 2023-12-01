#ifndef HH_LATTICE_READER_2D
#define HH_LATTICE_READER_2D

#include "latticeReader.hpp"
#include "lattice2D.hpp"
#include <filesystem>
#include <iostream>

/**
 * Implementation of the abstract LatticeReader class
 * used to read the sets of matrices that describe a 2D lattice and
 * the initialization values
*/
class LatticeReader2D : public LatticeReader
{
public:
    using LatticeGrid2D = typename Lattice2D::LatticeGrid2D;
    /**
     * Constructs the object
     * @param input_dir_path the path to the input directory
    */
    LatticeReader2D(const std::string input_dir_path_);
    virtual ~LatticeReader2D() = default;

    bool read_lattice_structure();
    bool read_lattice_input_rho();
    bool read_lattice_input_velocities();

};



#endif // HH_LATTICE_READER_2D