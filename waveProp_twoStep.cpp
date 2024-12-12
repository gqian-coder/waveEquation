#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include <dirent.h>
#include <random>
#include <algorithm>  // For std::min_element and std::max_element
#include <string.h>

#include "adios2.h"
#include "mgard/compress_x.hpp"
#include "waveEquation.hpp"
#include <zstd.h>
#include "waveInit.hpp"
#include "obstacle.hpp"

template <typename T>
void error_calc(std::vector<T> var_in, T *var_out, double tolerance)
{
    T diff, rmse=0.0;
    T minv = std::numeric_limits<T>::infinity();
    T maxv = -std::numeric_limits<T>::infinity();
    for (size_t i=0; i<var_in.size(); i++) {
        minv       = std::min(minv, var_in.data()[i]);
        maxv       = std::max(maxv, var_in.data()[i]);
        diff    = std::abs(var_in.data()[i] - var_out[i]);
        rmse   += diff*diff;
    }
    rmse      = std::sqrt(rmse / (T)var_in.size());

    //std::cout << "Error print out: \n";
    std::cout << "requested error: " << tolerance << ", ";
    std::cout << "L2: " << rmse << "\n"; 

}


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, np_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np_size);

    int cnt_argv = 1;
   
    // init_fun: 0 --> read from a data file
    //           1, 2, 3, 4, 5 --> sinusoidal, rain drop, multiple rains, square, velocity 
    // cfd_cond: 1, 2, 3 --> Dirichlet, Mur boundary condition, Nueman
    int init_fun   = std::stoi(argv[cnt_argv++]);
    int obstacle_t = std::stoi(argv[cnt_argv++]); 
    bool iterative_sim;
    std::istringstream(argv[cnt_argv++]) >> std::boolalpha >> iterative_sim;
    int cfd_cond   = std::stoi(argv[cnt_argv++]); 
    // simulation space 
    double Dx      = std::stof(argv[cnt_argv++]);
    double Dy      = std::stof(argv[cnt_argv++]);
    // simulation spatial resolution
    double dh      = std::stof(argv[cnt_argv++]);
    // simulation temporal resolution
    double dt      = std::stof(argv[cnt_argv++]);
    // number of simulation frames
    float T        = std::stof(argv[cnt_argv++]);
    // wave speed
    double C       = std::stof(argv[cnt_argv++]);
    // dissipation rate
    double gamma   = std::stof(argv[cnt_argv++]); 
    float temp     = T / dt;
    size_t nframes = (size_t)(temp);
    // relative tolerance. tol = 0 means to write out uncompressed data after each timestep
    double tol     = std::stof(argv[cnt_argv++]);
    // error distribution: assuming tol is for compressing u_curr, and tol * tol_ratio is for u_next
    double tol_ratio = std::stof(argv[cnt_argv++]);
    // error type
    std::string eb_type(argv[cnt_argv++]);
    // write out filename
    std::string fname_wt(argv[cnt_argv++]);
    // readin filename
    std::string fname(argv[cnt_argv++]);
    // initial timestep to read from file
    size_t init_ts = std::stoi(argv[cnt_argv++]);
    
    size_t Nx = (size_t)std::ceil((double)Dx / dh);
    size_t Ny = (size_t)std::ceil((double)Dy / dh);

    std::cout << "simulating a domain of [" << Dx << "/" << Nx << ", " << Dy << "/" << Ny << "], at a spacing of " << dh << "\n";
    std::cout << "simulating " << nframes << " steps at a resolution of " << dt << "\n";
    std::cout << "initial function: ";
    if (init_fun==0) std::cout << "read from previous simulation steps in the file: " << fname << " at step " << init_ts << "\n";
    else if (init_fun==1) std::cout << "sinusoidal \n";
    else if (init_fun==2) std::cout << "random raindrop\n";
    else if (init_fun==3) std::cout << "multiple rain drop\n";
    else if (init_fun==4) std::cout << "solid square\n";
    else if (init_fun==5) std::cout << "velocity fun in Gaussian\n";
    else if (init_fun==6) std::cout << "gaussian sinusoid waves\n";
    else if (init_fun==7) std::cout << "plane waves starting from x=0\n";
    else std::cout << "exponential wave velocity\n";
    std::cout << "simulating boundary condiction: ";
    if (cfd_cond==1) std::cout << "Dirichlet\n";
    else if (cfd_cond==2) std::cout << "Mur\n";
    else if (cfd_cond==3) std::cout << "Nueman\n"; 
    else if (cfd_cond==4) std::cout << "Periodic\n";
    std::cout << "dissipation rate = " << gamma << "\n";
    std::cout << "wave speed = " << C << "\n";
    std::cout << "tolerance = {" << tol << ", " << tol * tol_ratio << "} for u_n and u_dt \n";

    adios2::ADIOS ad(MPI_COMM_WORLD);
    adios2::IO writer_io = ad.DeclareIO("Output");

    adios2::Engine writer = writer_io.Open(fname_wt, adios2::Mode::Write);
    adios2::Variable<double> variable_u, variable_u_dt, variable_cr;
    variable_u  = writer_io.DefineVariable<double>("u_data" , adios2::Dims{Nx, Ny}, adios2::Dims{0,0}, adios2::Dims{Nx, Ny});
    std::vector<double> u_dt;
    if (tol>0) {
        // compression ratio of u_n and u_n-u_np1
        variable_cr = writer_io.DefineVariable<double>("u_CR"   , adios2::Dims{2}, adios2::Dims{0}, adios2::Dims{2});
        variable_u_dt = writer_io.DefineVariable<double>("u_dt", adios2::Dims{Nx, Ny}, adios2::Dims{0,0}, adios2::Dims{Nx, Ny});
        u_dt.resize(Nx*Ny);
    }
    /*
    adios2::Operator op = ad.DefineOperator("mgard", "mgard");
    if (tol>0){
        std::cout << "save the compressed data per timestep w/ eb = " << tol << "\n";
        variable_wt.AddOperation(op, {{"accuracy", std::to_string(tol)}, {"mode", "ABS"}}); 
    }
    */
    float max_intensity = 10.0;
    float freq = 2.0 * M_PI / T;
    float drop_probability = 1;
    std::vector<double> gauss_template;
    double NDx = 12, NDy=24;
    size_t n_drops = 5;
    double s1 = 1.0, s2 = 0.0; 

    // compression parameters
    mgard_x::Config config;
    config.lossless = mgard_x::lossless_type::Huffman_Zstd;
    //config.dev_type = mgard_x::device_type::SERIAL;
    config.dev_type = mgard_x::device_type::CUDA;
    //config.dev_id   = 1;
    std::vector<mgard_x::SIZE> shape{Nx, Ny};
    size_t compressed_size_u = 0, compressed_size_dt = 0;
    double data_bytes = (double) Nx * Ny * sizeof(double);
    std::vector<double> compression_ratio(2);
    WaveEquation <double> waveSim(Nx-1, Ny-1, 0, dt, dh, C, gamma, cfd_cond); 
    double *obstacle_m = NULL;
    if (obstacle_t) {
        switch (obstacle_t) {
            case 1: {
                size_t n_disks    = 6;
                size_t max_radius = size_t(0.02 * Nx);  
                obstacle_m        = emulate_N_disk<double>(Nx, Ny, n_disks, max_radius);
                break;
            }
            case 2: {
                obstacle_m = emulate_two_slits<double>(Nx, Ny);
                break;
            }
        }
    }

    if ((init_fun==2) || (init_fun==3)) {
        // Width of the Gaussian profile for each initial drop.
        double drop_width = 6;
        // Size of the Gaussian template each drop is based on.
        NDx = (size_t) std::ceil(drop_width / dh);
        NDy = (size_t) std::ceil(drop_width / dh);
        size_t cx = (size_t)(NDx/2);
        size_t cy = (size_t)(NDy/2);
        gauss_template.resize(NDx*NDy);
        std::vector <double> px(NDx), py(NDy);
        for (size_t r=0; r<NDx; r++) {
            px[r] = ((double)r - cx)/drop_width;
        }
        for (size_t c=0; c<NDy; c++) {
            py[c] = (double(c)-cy)/drop_width;
        }
        std::cout << "rainDrop region: " << NDx << " x " << NDy << " pixels\n";
        for (size_t r=0; r<NDx; r++) {
            for (size_t c=0; c<NDy; c++) {
                gauss_template[r*NDy+c] = max_intensity * exp(-(px[r]*px[r] + py[c]*py[c]));
            }
        }
    }
    switch (init_fun) {
        case 0: {
            adios2::IO reader_io = ad.DeclareIO("Input");
            reader_io.SetEngine("BP");
            adios2::Engine reader = reader_io.Open(fname, adios2::Mode::ReadRandomAccess);
            adios2::Variable<double> variable;
            if ((tol == 0) && (init_ts>0)) { // checkpoint restart 
                variable = reader_io.InquireVariable<double>("u_dt"); 
            } else {
                variable = reader_io.InquireVariable<double>("u_data");
            }
            std::cout << "total number of steps: " << variable.Steps() << ", read from " << init_ts << " timestep \n";
            std::vector<double> init_vpre(Nx*Ny);
            std::vector<double> init_v(Nx*Ny); 
            variable.SetSelection({adios2::Dims{0,0}, adios2::Dims{Nx, Ny}});
            variable.SetStepSelection({init_ts, 1}); 
            reader.Get(variable, init_vpre.data());
            reader.PerformGets();
            std::cout << "u_prev: min/max = " << *std::min_element(init_vpre.begin(), init_vpre.end()) << ", " << *std::max_element(init_vpre.begin(), init_vpre.end()) << "\n";
            std::cout << "begin to load the current step...\n";
            if ((init_ts>0)  && (tol==0)) {  
                variable = reader_io.InquireVariable<double>("u_data");
                variable.SetStepSelection({init_ts, 1});
            } else {
                variable.SetStepSelection({init_ts+1, 1});
            }
            reader.Get(variable, init_v.data());
            reader.PerformGets();
            std::cout <<  "u_curr: min/max = " << *std::min_element(init_v.begin(), init_v.end()) << ", " << *std::max_element(init_v.begin(), init_v.end()) << "\n";
            if ((init_ts>0) && (tol==0)) {  // change u_n - u_nm1 to u_mn1
                std::transform(init_v.begin(), init_v.end(), init_vpre.begin(), init_vpre.begin(), std::minus<double>());
            }
            waveSim.init_u_n(init_vpre.data());
            waveSim.init_u_np1(init_v.data());
            reader.Close();
            init_v.clear();
            init_vpre.clear(); 
            if (tol==0) {
                writer.BeginStep();
                writer.Put<double>(variable_u, waveSim.u_n.data(), adios2::Mode::Sync);
                writer.PerformPuts();
                writer.EndStep();
                writer.BeginStep();
                writer.Put<double>(variable_u, waveSim.u_np1.data(), adios2::Mode::Sync);
                writer.PerformPuts();
                writer.EndStep();
            }
            break;
        }
        case 1:
            break;
        case 2: { 
            fun_rainDrop<double>(waveSim.u_np1.data(), Nx, Ny, NDx, NDy, gauss_template.data(), drop_probability);
            break;
        }
        case 3: {
            fun_MultiRainDrop<double>(waveSim.u_np1.data(), Nx, Ny, NDx, NDy, gauss_template.data(), drop_probability, n_drops);
            break;
        }
        case 4: {
            fun_square<double>(waveSim.u_np1.data(), Nx, Ny, NDx, NDy, max_intensity);
            break;
        }
        case 5: {
            double sigma = 0.1;
            fun_velocity<double>(waveSim.u_np1.data(), Nx, Ny, dh, dh, max_intensity, sigma);
            break;
        }
        case 6: {
            int n_waves = 4;
            std::vector<double> sigma = {5, 2, 7.5, 10};
            std::vector<double> intensity = {25, 30, 45, 15};
            std::vector<double> freq_x = {freq*60, freq*85, freq*95, freq*50};
            std::vector<double> freq_y = {freq*70, freq*85, freq*105, freq*68};
            fun_gaussian_wave(waveSim.u_np1.data(), Nx, Ny, intensity, sigma, freq_x, freq_y, n_waves);
            break;
        }
        case 7: {
            size_t n_waves = size_t(0.02 * (double)Ny);
            fun_plane_waves<double>(waveSim.u_np1.data(), Nx, Ny, max_intensity, freq, n_waves);
            break;
        }
        default:
            break; 
    }
    if (obstacle_t) {
            // apply Dirichlet boundary condition to the scattering on obstacles
            std::transform(waveSim.u_np1.begin(), waveSim.u_np1.end(), obstacle_m, waveSim.u_np1.begin(), std::multiplies<double>());
    }

    drop_probability  = 0.01; 
    size_t iter_frame = 0;
    while (iter_frame<nframes) {
        std::cout << iter_frame << "/" << nframes << "\n";
        writer.BeginStep();
        // wave update 
        waveSim.update_2d(static_cast<bool>(iter_frame)); 
        if (obstacle_t) {
            // apply Dirichlet boundary condition to the scattering on obstacles
            std::transform(waveSim.u_np1.begin(), waveSim.u_np1.end(), obstacle_m, waveSim.u_np1.begin(), std::multiplies<double>());
        }
        if (tol>0) {
            void *compressed_array_cpu  = NULL;
            void *compressed_array2_cpu = NULL;
            if (strcmp(eb_type.c_str(), "ABS")==0) {
                mgard_x::compress(2, mgard_x::data_type::Double, shape, tol, s1,
                    mgard_x::error_bound_type::ABS, waveSim.u_np1.data(),
                    compressed_array_cpu, compressed_size_u, config, false);
                
                std::transform(waveSim.u_np1.begin(), waveSim.u_np1.end(), waveSim.u_n.begin(), u_dt.begin(), std::minus<double>());
                mgard_x::compress(2, mgard_x::data_type::Double, shape, tol*tol_ratio, s2,
                    mgard_x::error_bound_type::ABS, u_dt.data(),
                    compressed_array2_cpu, compressed_size_dt, config, false);
            } else {
                mgard_x::compress(2, mgard_x::data_type::Double, shape, tol, s1,
                    mgard_x::error_bound_type::REL, waveSim.u_np1.data(),
                    compressed_array_cpu, compressed_size_u, config, false);

                std::transform(waveSim.u_np1.begin(), waveSim.u_np1.end(), waveSim.u_n.begin(), u_dt.begin(), std::minus<double>());
                mgard_x::compress(2, mgard_x::data_type::Double, shape, tol*tol_ratio, s2,
                    mgard_x::error_bound_type::REL, u_dt.data(),
                    compressed_array2_cpu, compressed_size_dt, config, false);
            }
            compression_ratio[0] = data_bytes / (double) compressed_size_u;
            compression_ratio[1] = data_bytes / (double) compressed_size_dt;
            std::cout << "compression ratios = " << compression_ratio[0] << " and " << compression_ratio[1] << "\n";
            void *decompressed_array_cpu = NULL;
            void *decompressed_array2_cpu = NULL;
            mgard_x::decompress(compressed_array_cpu, compressed_size_u,
                decompressed_array_cpu, config, false);
            mgard_x::decompress(compressed_array2_cpu, compressed_size_dt,
                decompressed_array2_cpu, config, false);
            
            //std::cout << ((double *)decompressed_array_cpu)[1000] << ", " << ((double *)decompressed_array2_cpu)[1000] << "\n";
            error_calc(waveSim.u_np1, (double *)decompressed_array_cpu, tol); 
            error_calc(u_dt, (double *)decompressed_array2_cpu, tol*tol_ratio);
            writer.Put<double>(variable_u   , (double *)decompressed_array_cpu , adios2::Mode::Sync);
            writer.Put<double>(variable_u_dt, (double *)decompressed_array2_cpu, adios2::Mode::Sync);
            writer.Put<double>(variable_cr, compression_ratio.data(), adios2::Mode::Sync);
            if (iterative_sim) {
                std::cout << "iterative update solution over time step\n";
                waveSim.init_u_np1((double *)decompressed_array_cpu);
            }
        } else {
            writer.Put<double>(variable_u, waveSim.u_np1.data(), adios2::Mode::Sync);
        }
        iter_frame ++;
        writer.PerformPuts();
        writer.EndStep(); 
    }
    writer.Close();
    
    delete [] obstacle_m;
    MPI_Finalize();
    return 0;
}













