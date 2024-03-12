# Set the source directory and gather source files
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Create llalbm target
add_library(LLALBM INTERFACE)

# Find dependencies: MPI, OpenMP, and Eigen
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE HINTS $ENV{mkEigenPrefix})

# Link libraries to the llalbm target
target_include_directories(LLALBM INTERFACE 
	"${CMAKE_CURRENT_SOURCE_DIR}/include"
	"${CMAKE_CURRENT_SOURCE_DIR}/src"
	"${EIGEN3_INCLUDE_DIR}")

message("LLALBM INCLUDE DIRECTORY : ${CMAKE_CURRENT_SOURCE_DIR}/include")
