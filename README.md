# LLALBM
In this repository we present __LLALBM__, a templated and header only C++ library implementing the Lattice Boltzmann Method for computational fluid dynamics.

## Cloning, Compiling and Executing
### Before Cloning
Before working with LLALBM, please make sure you have all the libraries we depend on correctly installed in your system. You will need
- __Eigen__
- __OpenMP__
- __OpenACC__
- __MPI__
- __CMAKE__

We used __mk modules__ by the __MOX__ laboratory of __Politecnico di Milano__. You can find them [here](https://github.com/pcafrica/mk).
### Cloning
In order to work with LLALBM, start by cloning the repository.
```shell
git clone https://github.com/AMSC22-23/LatticeBoltzmannMethod-Fonnesu-Grassi-Guffanti.git
```
or
```shell
git clone --recurse-submodules git@github.com:AMSC22-23/LatticeBoltzmannMethod-Fonnesu-Grassi-Guffanti.git
```
Once you have cloned the repository, you are ready to compile you first examples!

### Compiling
Compilation is supported by [CMake](https://cmake.org/).
In the [CMakeLists.txt](/CMakeLists.txt) in the root of the repository you can configure the compilation of LLALBM by 
choosing whether to enable compilation of tests by acting on the `should_Test` option. Code examples are compiled by default and
are accessible in the [examples](/examples/) directory.

Once you have chosen how to compile LLALBM, you are ready to actually build it. To do so,
create a `build/` folder in the root of the repository and move into it
```shell
mkdir build
cd build/
```
Remember now to load all the modules you need, and from 
inside the `build/` folder simply call CMake to produce Makefiles, and then compile the project with `make`.
```shell
cmake ..
make
```
Please notice that based on your LLALBM configuration compilation may take different amounts of time.
### Executing LLALBM examples
Once compilation is finished, you can easily run examples. If you have compiled the examples from `build/`, you
will find examples in the `build/examples/` directory. If you have compiled the tests too, you find executables in the
`build/test/` folder. 

## Example Gallery.
LLALBM comes packaged with a number of example that show its capabilities, and whose source code can be easily accessed from
the [examples](/examples/) directory. C++ files containing the example code have the same 
name of the example you find right here.
### [1_LidDrivenCavity](/examples/1_LidDrivenCavity.cpp) and [1_LidDrivenCavity129](/examples/1_LidDrivenCavity129.cpp)

### [2_Parabolic](/examples/2_Parabolic.cpp)

### [3_Vortex](/examples/3_Vortex.cpp)

### [4_ComplexGeneration](/examples/4_ComplexGeneration.cpp)

### [5_RandomGeneration](/examples/5_RandomGeneration.cpp)

### [6_FlowAnalysis](/examples/6_FlowAnalysis.cpp)

### [7_FlatBottomFoil](/examples/7_FlatBottom)

### [8_1NACAProfileFlight](/examples/8_1NACAProfileFlight.cpp)

### [8_2NACAProfileTakeOff](/examples/8_2_NACAProfileTakeOff.cpp)

### [9_Conduct](/examples/9_Conduct.cpp)

## Visualising the results yourself

Implemented by
- [Lorenzo Fonnesu](mailto:lorenzo.fonnesu@mail.polimi.it)
- [Andrea Grassi](mailto:andrea10.grassi@mail.polimi.it)
- [Luca Guffanti](mailto:luca2.guffanti@mail.polimi.it)
