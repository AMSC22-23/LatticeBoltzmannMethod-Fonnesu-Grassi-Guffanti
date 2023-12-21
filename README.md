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
### Computing the lattice
Lattices describe computational domains and are stored in `.mtx` files. We provide a few example lattices that can be used as well as a few `.png` files from which lattices can be generated. To generate a lattice you'll need to call the `translate.py` python script in the following way.
```shell
python ./scripts/translate.py path_to_image lattice path_to_output_dir
```
### Executing the code
To run the code you'll simply need to call, within `build`, `./lattice_boltzmann-seq` or `./lattice_boltzmann-omp` with the adeguate set of parameters. Both executables accept the same parameters.
```shell
./lattice_boltzmann-(seq|omp) dimensions input_dir collision_model [-r reynolds] [-f output freq]
```
- `dimensions` 2 or 3, describes the dimensions of the simulation. Until now we only support 2D simulations.
- `input_dir` path to the input directory containing the lattice.mtx file.
- `collision model` only **BGK** or **TRT** are available: this describes the model used to execute collision between fluid nodes.
- `reynolds` the Reynolds number of the simulation. Beware that high Reynolds number may cause numerical instability.
- `output freq` frequency of the output file creation in the case of a simulation (strong scaling does not produce files).
### Visualizing the results
You have two python scripts at your disposal to visualize the results, they are in the `scripts/` directory.

# Credits
Implemented by
- [Lorenzo Fonnesu](mailto:lorenzo.fonnesu@mail.polimi.it)
- [Andrea Grassi](mailto:andrea10.grassi@mail.polimi.it)
- [Luca Guffanti](mailto:luca2.guffanti@mail.polimi.it)

as a project for the Advanced Methods for Scientific Computing course at Politecnico di Milano
# Bibliography and Used Resources
[1] U. Ghia, K. N. Ghia, C. T. Shin, High-Re solutions for incompressible flow using Navier-Stokes equations
and multigrid method.
