# LLALBM, A Lightweight, Lightspeed and Adaptable LBM implementation.
In this repository we present __LLALBM__, a templated and header only C++ library implementing the Lattice Boltzmann Method for computational fluid dynamics

<blockquote style="color: red;">
This README provides a description of how to clone, compile, execute and visualize simulations with LLALBM. For a comprehensive description of the
architecture of the project, both from the numerical and software engineering viewpoints, please refer to the report you find in the root of the repository.
</blockquote> 



## Cloning, Compiling and Executing
### Before Cloning
Before working with LLALBM, please make sure you have all the libraries we depend on correctly installed in your system. You will need
- __Eigen__
- __OpenMP__
- __OpenACC__
- __MPI__
- __CMAKE__

We used __mk modules__ by the __MOX__ laboratory of __Politecnico di Milano__. You can find them [here](https://github.com/pcafrica/mk).

<blockquote style="color: red;">
Modules are preferred over git submodules thanks to their flexibility, and actual use in high performance systems.
</blockquote>  

### Cloning
In order to work with LLALBM, start by cloning the repository.
```shell
git clone https://github.com/AMSC22-23/LatticeBoltzmannMethod-Fonnesu-Grassi-Guffanti.git
```
or
```shell
git clone git@github.com:AMSC22-23/LatticeBoltzmannMethod-Fonnesu-Grassi-Guffanti.git
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
The Lid Driven Cavity example is the basic computational fluid dynamics benchmark. We impose a horinzontal flow modeling
a moving wall on the top boundary of the domain. The two files represent simulations with two different sizes. This first example is
used to show how easy it is to generate a computational domain with the built-in lattice generation system.
![LidDrivenGif](/examples/videos/example_1.gif)
The Lid driven simulation complies with the suggested literature, as shown by the following plots.
<img src="/examples/data/images/velocities-1.png" alt="LidDrivenVx" width="500">
<img src="/examples/data/images/v-velocities-1.png" alt="LidDrivenVy" width="500">
### [2_Parabolic](/examples/2_Parabolic.cpp)
The second example showcases a flow near the turbulent regime in a channel with a parabolic inlet flow profile in the horizontal direction.
Such simulation is useful to analyze how pivotal the parameter choice is: the simulation has in fact some artifacts, most probably due to the
introduction of numerical noise.
![ParabolicGif](/examples/videos/example_2.gif)
### [3_Vortex](/examples/3_Vortex.cpp)
The third example is used to show how to introduce an obstacle, in this case a cylinder, in a computational domain. The cylinder is placed
asymmetrically in the domain in order to stimulate, especially in the case of a high Reynolds number, a vortex street behind the cylinder. Two simulations
are put in place, the first with a viscous fluid (altough a hint of vortices is still visible), and the second with a less viscous fluid.
![VortexVisc](/examples/videos/example_3_LOWRE_1.gif)
![VortexNotVisc](/examples/videos/example_3_HIGHRE.gif)

### [4_ComplexGeometry](/examples/4_ComplexGeometry.cpp) 
LLALBM in-code domain generation is very well suited for generating lattices with a big number of obstacles. Up to now LLALBM supports natively
the generation of Rectangles, Squares, Spheres. Nonetheless, functionalities exposed by the lattice reading infrastructure allow the reading
of obstacle information from an input stream, greatly extending the range of constructable domains. (Lattice reading to introduce obstacle is 
showcased in the following examples.)
![Complexgif](/examples/videos/example_4.gif)

### [5_RandomGeneration](/examples/5_RandomGeneration.cpp)
As a way of trying to simulate how a flow would behave in a porous media, or in a solid media where nucleation phenomena occurred, LLALBM supports 
a random domain generation. Thus, this example shows how a random domain can be built by using this functionality, with the capability of creaing interesting flows.
Random number generation employs the random number generation infrastructure put in place by the C++ standard library
![Randomgif](/examples/videos/example_5_Random_1.gif)

### [6_FlowAnalysis](/examples/6_FlowAnalysis.cpp)
LLALBM also features a flow analysis subsystem (centered around the `FlowAnalyzer<dim>` object) that evaluates lift and drag forces on a set 
of obstacle nodes chosen by the user. The lift analysis can be executed together with the simulation or after the simulation. This example can be configured
to execute in both ways. The domain contains a cylinder with a rectangular tail placed in the center to avoid the generation of the vortex street after the obstacle.
Results of the simulation are stored `.csv` files and can be rendered into plots with the provided python scripts.
![FlowAnalysisgif](/examples/videos/example_6.gif)
<img src="/examples/data/images/6_Drag.png" alt="drag6" width="500">
<img src="/examples/data/images/6_Lift.png" alt="lift6" width="500">
<img src="/examples/data/images/6_Lift_Drag-1.png" alt="liftdrag6" width="500">
### [7_FlatBottomFoil](/examples/7_FlatBottom)
After considering all the previous results, we can now simulate a more complex flow: the one around an aerofoil. This example implements
a simulation of a Flat Bottom aerofoil.
The simulation is used to produce an animation and a lift and drag analysis is put in place to understand whether the use of a Flat Bottom aerofoil is beneficial.

![FlatBottomAeroGifLaminar](/examples/videos/example_7.gif)
<img src="/examples/data/images/7_Lift_new.png" alt="drag7" width="500">
<img src="/examples/data/images/7_drag_new.png" alt="lift7" width="500">
<img src="/examples/data/images/7_Lift_Drag_new-1.png" alt="liftdrag7" width="500">

### [8_1NACAProfileFlight](/examples/8_1NACAProfileFlight.cpp)
In this example a NACA aerofoil is explored in normal flight conditions in two different viscosity conditions. A uniform horizontal 
velocity is enforced on the inlet nodes that present on the left of the domain. Additionally, a lift analysis is carried out to evaluate the behaviour of 
the aerofoil in the flow.

The above simulation is executed in a fluid with high viscosity, while the one below happens in a less viscous fluid. 
![NacaFlightLaminar](/examples/videos/example_8_1_LOWRE.gif)
![NacaFlightTurbulent](/examples/videos/example_8_1_HIGHRE.gif)

And the flow analysis is presented here: on the left the forces for the viscous fluid, and on the right the forces for the smooth fluid.
<img src="/examples/data/images/8_1_Lift_Drag_LOWRE_new-1.png" alt="draglift8" width="500">
<img src="/examples/data/images/8_1_lift_drag_HIGHRE_new-1.png" alt="draglift82" width="500">


### [8_2NACAProfileTakeOff](/examples/8_2_NACAProfileTakeOff.cpp)
In this example, which completes the previous one, a NACA aerofoil is simulated during takeoff. To address the different takeoff conditions, the mesh is modified and rotated
to simulate a variation in angle of attack. In a first phase a uniform horizontal velocity is imposed on the inlet nodes to simulate taxiing, and then vertical displacement is simulated by introducing vertically flowing fluid from the lower part of the domain boundary. 
Again, a lift and drag analysis is performed during the simulations, which are constructed in a laminar regime and in a turbulent regime.

The above simulation is executed in a fluid with high viscosity, while the one below is in a less viscous fluid.
![NacaTakeOffLaminar](/examples/videos/example_8_2_LOWRE.gif)
![NacaTakeOffTurbulent](/examples/videos/example_8_2_HIGHRE.gif)

The results of the analysis are presented here, with the laminar regime on the left and the turbulent regime on the right.
<img src="/examples/data/images/8_2_Lift_Drag_LOWRE_new-1.png" alt="draglift812" width="500">
<img src="/examples/data/images/8_2_lift_drag_HIGHRE_new-1.png" alt="draglift822" width="500">


### [9_Conduct](/examples/9_Conduct.cpp)
Using again the obstacle generation, LLALBM can be used to simulate fluid motion in conducts. Conducts can be read from file or directly generated 
in the code. The example is shown below
![ConductGif1](/examples/videos/example_9_1.gif)
![ConductGif](/examples/videos/example_9_norm2.gif)

## Visualising the results yourself
In order to build the simulations and plots that are presented in the previous sections, you will need to use some __Python__ scripts!
However, our scripts depend on a few libraries you will need to install. The approach we suggest is constructing a virtual environment in which 
you can freely download packages without the risk of collisions with already present packages.
### Installing Python environment and packages
First, make sure you have __Python>=3.3__, as virtual environments have been introduced by that release.

Then, build and a virtual environment with
```shell
python -m venv path_to_venv_directory
```
and activate it with
```shell
source path_to_venv_directory/bin/activate
```
To install packages you will just need pip, which is distributed with Python.
First, upgrade pip with 
```shell
python -m pip install --upgrade pip 
```
and you can now download packages by running
```shell
pip install -r requirements.txt 
```
You should now be up and running!
### Python scripts
In this repository, in the [scripts](/scripts) folder, you find a few scripts that can be used
to generate lattices and obstacles, render animations, and plot the flow analysis.

-`translate.py path_to_image lattice|rho path_to_output_dir outputfilename` can be used to generate lattices.


-`produce_video.py path_to_input_dir u (norm)` can be used to generate animations starting from velocity fields.


-`plot_obstacle_layout_2d.py path_to_file` is used to plot an obstacle layout from a .csv file generated from FlowAnalyzer. 


-`plot_lift_drag.py path_to_dir` is used to plot single nodal contributions as well as the global flow and lift. 


-`plot_data_x.py file_path_x file_path_y line` is used to plot velocities on a horizontal line.


-`plot_data_y.py file_path_x file_path_y line` is used to plot velocities on a vertical line.


-`build_obstacle_file.py input_file output_file` is used to traduce a .png image representing an obstacle into a .txt
file comprehensible by LLALBM.


## Acknowledgements
Implemented by
- [Lorenzo Fonnesu](mailto:lorenzo.fonnesu@mail.polimi.it)
- [Andrea Grassi](mailto:andrea10.grassi@mail.polimi.it)
- [Luca Guffanti](mailto:luca2.guffanti@mail.polimi.it)
