#include "lbm.hpp"

lbm::lbm(const std::size_t& D,
    const double reynolds_, 
    const std::string& input_dir_path, 
    const std::string& collision_model,
    const std::size_t& frequency_ = 1):
    frequency (frequency_),
    re (reynolds_)
    {
    std::string output_dir_path="../results";
    lattice_ptr=nullptr;

    constexpr double one_third = 1.0/3.0;
    const double tau = 0.5 + (0.2 * 100.0 /reynolds_)/(one_third);
    const double t_const = 1.0 / tau;
    const double t_conj = 1.0 - t_const;

    std::cout << "tau   : " << tau << std::endl;
    std::cout << "t/tau : " << delta_t/tau << std::endl;
    std::cout << "t_const : " << t_const << std::endl;
    std::cout << "t_conj : " << t_conj << std::endl;

    if(collision_model == "BGK")
    {
        collision_ptr = std::make_shared<BGK>(t_conj, t_const);
    }else if(collision_model == "TRT")
    {
        collision_ptr = std::make_shared<TRT>(tau);
    }else if(collision_model == "MRT")
    {
        collision_ptr = std::make_shared<MRT>();
    }else
    {
        std::cout << collision_model << " not yet implemented" << std::endl;
    }

    switch (D)
    {
        case 2:
            velocity_set.initialize(2, 9);
            lattice_ptr = std::make_unique<Lattice2D>(input_dir_path, output_dir_path, velocity_set, collision_ptr, re);
            break;
        case 3:
            lattice_ptr=nullptr;
            break;
        default:
            std::cout << "[!ERROR!] only 2D and 3D are supported" << std::endl;
            break;
    }
}

void lbm::compute(const double time)
{
    double size = 500.0;
    double dt = re * (0.2 * size /re)/(size * size);
    const std::size_t n_iter = floor(time/dt);

    std::cout << "Time  : " << time << std::endl;
    std::cout << "dt    : " << dt << std::endl;
    std::cout << "iter  : " << n_iter << std::endl;

    std::size_t output_counter = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i <= n_iter; i++)
    {
        lattice_ptr->set_inlets(i);
        // REINITIALIZATION OF INLET FIELDS
        lattice_ptr->perform_simulation_step();
        if(i % frequency == 0)
        {
            
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            std::cout << duration.count() << " ms" << std::endl;
            // print the output only if the iteration is a multiple of the frequency. 
            lattice_ptr->save_output_data(output_counter);
            output_counter++;
            start_time = std::chrono::high_resolution_clock::now();
        }
    }    
}

void lbm::perform_strong_scaling_test()
{
    #ifdef _OPENMP
        double dt = re * (0.2 * 100 /re)/(100.0 * 100.0);
        const std::size_t n_iter = floor(5.0/dt);
        std::vector<std::tuple<std::size_t, double>> results_table;

        std::size_t max_num_threads = omp_get_max_threads();
        std::cout << "=====STRONG SCALING TEST=====" << std::endl;
        std::cout << "Available OMP threads: " << max_num_threads << std::endl;

        for (std::size_t num_threads = 1; num_threads <= max_num_threads; num_threads++)
        {
            std::cout << "Testing with " << num_threads << " thread(s): ";
            lattice_ptr->set_omp_num_threads() = num_threads;
            auto start_time = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i <= n_iter; i++)
            {
                {
                    // REINITIALIZATION OF INLET FIELDS
                    lattice_ptr->set_inlets(i);

                    lattice_ptr->perform_simulation_step();                    
                }    
            }
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            std::cout << duration.count()/1000.0 << " ms" << std::endl;
            
            results_table.push_back({num_threads, duration.count()/1000.0});
        }
        std::ofstream out_file("strong_scalability.txt");
        for (const auto& el : results_table)
        {
            out_file << std::get<0>(el) << " " << std::get<1>(el) << std::endl;
        }
        out_file.close();
    #else
        std::cout << "[!ERROR!] Strong scalability test is unavailable: compile with -fopenmp" << std::endl;
    #endif
}


void lbm::perform_weak_scaling_test(size){
    #ifdef _OPENMP

        std::vector<std::tuple<int,std::size_t, double>> results_table;
        std::size_t max_num_threads = omp_get_max_threads();
        
        std::cout << "=====WEAK SCALING TEST=====" << std::endl;
        std::cout << "Available OMP threads: " << max_num_threads << std::endl;
    
        for (std::size_t num_threads = 1; num_threads <= max_num_threads; num_threads++)
        {   
            double dt = re * (0.2 * size /re)/(size * size);
            const std::size_t n_iter = floor(5.0/dt);
            std::cout << "Testing with " << num_threads << " thread(s): ";

            
            lattice_ptr->set_omp_num_threads() = num_threads;
            auto start_time = std::chrono::high_resolution_clock::now();
            for (std::size_t i = 0; i <= n_iter; i++)
            {
                {
                    // REINITIALIZATION OF INLET FIELDS
                    lattice_ptr->set_inlets(i);

                    lattice_ptr->perform_simulation_step();                    
                }    
            }
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            std::cout << duration.count()/1000.0 << " ms" << std::endl;
            
            results_table.push_back({size ,num_threads, duration.count()/1000.0});
        }
    

        std::ofstream out_file("weak_scalability" +size+".txt");
        for (const auto& el : results_table)
        {
            out_file << std::get<0>(el) << " " << std::get<1>(el) << " " << std::get<2>(el) << std::endl;
        }
        out_file.close();
        
    #else
       std::cout << "[!ERROR!] Weak scalability test is unavailable: compile with -fopenmp" << std::endl;
    #endif
}
