# Lattice Boltzmann method
The Lattice Boltzmann method is a mesoscopic method for computational fluid dynamics arising from the discretization of the Boltzmann equation: one of the most known results of statistical mechanics, formulated by Ludwig Boltzmann in 1872, that describes the statistical behaviour of a fluid in a non-equilibrium state.
## Scope of the project
We implemented a C++ & OpenMP version of the method, currently supporting the lid-driven cavity scenario: a benchmark problem for viscous incompressible fluid flow. To complement our implementation we performed weak and strong scaling tests to check the quality of the parallelization and a correctness test against the results from [1] considering the center of the computational domain.
# Structure of the repository
### Directory structure
### Use of submodules
Executing the code requires [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) that we included as a git submodule. We chose to include eigen as a git submodule to allow for ease of use of the codebase without the need of using mk-modules and to always have the updated version. However, you're free to use mk-modules! We use Eigen for its support of matrices (and in the near future for vectors and tensors, aiming also to port all native std::vectors/arrays to Eigen vectors).

# Cloning, Compiling and Executing
### Cloning


### Compiling
Compilation is supported by [CMake](https://cmake.org/). To compile the project, first create a `build` directory and move into it
```shell
mkdir build
cd build/
```
Then call CMake to build the makefile and subsequently compile the project.
```shell
cmake ..
make
```
You'll see all the compilation log, and two executables will be created: `lattice_boltzmann-seq` and `lattice_boltzmann-omp`. The code is basically the same, with `lattice_boltzmann-omp` that supports multithreading with OpenMP and allows execution of the weak scaling test. 
### Executing the code
To run the code you'll simply need to call, within `build`, `./lattice_boltzmann-seq` or `./lattice_boltzmann-omp`. 
### Visualizing the results

# Credits
Implemented by
- Lorenzo Fonnesu
- Andrea Grassi
- Luca Guffanti

as a project for the Advanced Methods for Scientific Computing course at Politecnico di Milano
# Bibliography and Used Resources
[1] U. Ghia, K. N. Ghia, C. T. Shin, High-Re solutions for incompressible flow using Navier-Stokes equations
and multigrid method.
