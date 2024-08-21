// =========== STL INCLUDES ===========
#include <vector>
#include <array>
#include <cassert>
#include <functional>
// ======================================

// =========== EIGEN INCLUDES ===========
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
// ======================================

// =========== CUDA INCLUDES ===========
#include <cuda_runtime.h> 
#include <device_launch_parameters.h>
// ======================================

__global__ void update_macro_kernel(const double* populations, const Point<2>* fluid_nodes, 
                                    double* global_rho, double* global_u, int num_fluid_nodes, int width, int height) {
    int fnode = blockIdx.x * blockDim.x + threadIdx.x;

    if (fnode < num_fluid_nodes) {

        int i = fluid_nodes[fnode].coords[0];
        int j = fluid_nodes[fnode].coords[1];

        int idx = i * height + j;

        double p0 = populations[idx * 9 + 0];
        double p1 = populations[idx * 9 + 1];
        double p2 = populations[idx * 9 + 2];
        double p3 = populations[idx * 9 + 3];
        double p4 = populations[idx * 9 + 4];
        double p5 = populations[idx * 9 + 5];
        double p6 = populations[idx * 9 + 6];
        double p7 = populations[idx * 9 + 7];
        double p8 = populations[idx * 9 + 8];

        double rho = p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8;
        double rhoinv = 1.0 / rho;

        double ux = rhoinv * (p1 + p5 + p8 - (p3 + p6 + p7));
        double uy = rhoinv * (p2 + p5 + p6 - (p4 + p7 + p8));

        global_rho[idx] = rho;  
        global_u[idx * 2 + 0] = ux;
        global_u[idx * 2 + 1] = uy;
    }
}

void update_macro(const Tensor<double, 3> &populations, const std::vector<Point<2>> &fluid_nodes, Tensor<double, 2> &global_rho, Tensor<double, 3> &global_u)
{    

    size_t num_fluid_nodes = fluid_nodes.size();
    size_t width = populations.dimension(0);
    size_t height = populations.dimension(1);
    
    // Aloocation of the memory on the GPU
    double *d_populations, *d_global_rho, *d_global_u;
    Point<2>* d_fluid_nodes;

    cudaMalloc(&d_populations, width * height * 9 * sizeof(double));
    cudaMalloc(&d_global_rho, width * height * sizeof(double));
    cudaMalloc(&d_global_u, width * height * 2 * sizeof(double));
    cudaMalloc(&d_fluid_nodes, num_fluid_nodes * sizeof(Point<2>));

    // Copy data from CPU to GPU
    cudaMemcpy(d_populations, populations.data(), width * height * 9 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fluid_nodes, fluid_nodes.data(), num_fluid_nodes * sizeof(Point<2>), cudaMemcpyHostToDevice);

    // Kernel size
    int blockSize = 256;
    int numBlocks = (num_fluid_nodes + blockSize - 1) / blockSize;

    // Kernel call
    update_macro_kernel<<<numBlocks, blockSize>>>(d_populations, d_fluid_nodes, d_global_rho, d_global_u, num_fluid_nodes, width, height);

    // Check for errors in the kernel
    cudaError_t kernelErr = cudaGetLastError();
    if (kernelErr != cudaSuccess) {
        std::cerr << "Errore nel kernel: " << cudaGetErrorString(kernelErr) << std::endl;
        cudaFree(d_populations);
        cudaFree(d_fluid_nodes);
        cudaFree(d_global_rho);
        cudaFree(d_global_u);
        return;
    }

    // Copy of results from GPU to CPU
    cudaMemcpy(global_rho.data(), d_global_rho, width * height * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(global_u.data(), d_global_u, width * height * 2 * sizeof(double), cudaMemcpyDeviceToHost);

    // Free memory
    cudaFree(d_populations);
    cudaFree(d_global_rho);
    cudaFree(d_global_u);
    cudaFree(d_fluid_nodes);
}