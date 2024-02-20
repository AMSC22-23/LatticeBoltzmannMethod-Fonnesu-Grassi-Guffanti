# Create llalbm target
add_library(LLALBM SHARED ${SOURCES})


# include OpenMP, MPI and Eigen
find_package(MPI REQUIRED)
cmake_policy(SET CMP0074 NEW)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE
		HINTS $ENV{mkEigenPrefix})


# Link libraries to the executable
target_link_libraries(LLALBM PUBLIC MPI::MPI_CXX Eigen3::Eigen OpenMP::OpenMP_CXX)

# Additional settings for MPI
target_include_directories(LLALBM PUBLIC ${MPI_C_INCLUDE_PATH})