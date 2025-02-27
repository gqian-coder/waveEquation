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

// potential energy
double calc_PE(double *u, double *pe, double dh, size_t Nx, size_t Ny)
{
    double PE = 0.0;
    double ux, uy;
    size_t r_curr, k;
    for (size_t r=1; r<Nx; r++) {
        r_curr = r * Ny;
        for (size_t c=1; c<Ny; c++) {
            k     = r_curr + c;
            ux    = (u[k] - u[k-Ny]) / dh;
            uy    = (u[k] - u[k-1])/dh;
            pe[k] = ux*ux + uy*uy;
            PE   += pe[k];
        }
    }
    return PE * 0.5;
}

// root-of-mean-square error
double calc_rmse(double *data_f, double *data_g, double *diff, size_t num_data)
{
    double rmse = 0.0;
    for (size_t i=0; i<num_data; i++) {
        diff[i] = data_f[i] - data_g[i];
        rmse += (diff[i]*diff[i]) / (double) num_data; 
    }
    return std::sqrt(rmse); 
}

int main(int argc, char **argv) {

    int cnt_argv = 1;
    std::string fname_f(argv[cnt_argv++]);
    // simulation spatial resolution
    double dh      = std::stof(argv[cnt_argv++]);
    // from the init_ts step in frame_f to compare against frame_g 
    size_t init_ts = std::stoi(argv[cnt_argv++]);    
    std::string fname_err(argv[cnt_argv++]);
    double tol_s0  = std::stof(argv[cnt_argv++]);
    double tol_s1  = std::stof(argv[cnt_argv++]);  

    std::cout << "original file: " << fname_f.c_str() << "\n";
    std::cout << "output file: " << fname_err.c_str() << "\n"; 
    std::cout <<  "dh = " << dh << ", init_ts = " << init_ts << ", tol(s0) = " << tol_s0 << ", tol(s1) = " << tol_s1 << "\n";

    adios2::ADIOS ad;
    adios2::IO reader_io_f = ad.DeclareIO("Original");
    reader_io_f.SetEngine("BP");
    adios2::Engine reader_f = reader_io_f.Open(fname_f, adios2::Mode::ReadRandomAccess);
    adios2::Variable<double> variable_f = reader_io_f.InquireVariable<double>("u_data");
    
    size_t Nx = variable_f.Shape()[0];
    size_t Ny = variable_f.Shape()[1];
    size_t num_data = Nx * Ny; 
    std::vector<double> var_f(num_data);
    // difference data
    std::vector<double> err_s1(num_data);
    std::vector<double> err_s0(num_data);
    std::vector<double> PE_s0(num_data);
    std::vector<double> PE_s1(num_data);

    // compression parameters
    mgard_x::Config config;
    config.lossless = mgard_x::lossless_type::Huffman_Zstd;
    //config.dev_type = mgard_x::device_type::SERIAL;
    config.dev_type = mgard_x::device_type::CUDA;
    std::vector<mgard_x::SIZE> shape{Nx, Ny};

    variable_f.SetStepSelection({init_ts, 1});
    reader_f.Get(variable_f, var_f);
    reader_f.PerformGets();
    
    double v_max = *std::max_element(var_f.begin(), var_f.end());
    double v_min = *std::min_element(var_f.begin(), var_f.end());
    std::cout << "data value ranges: {" << v_max << ", " << v_min << "}\n";

    void *compressed_s0   = NULL, *compressed_s1   = NULL;
    void *decompressed_s0 = NULL, *decompressed_s1 = NULL;
    size_t compressed_size_s0, compressed_size_s1;
    // s=0 compression
    mgard_x::compress(2, mgard_x::data_type::Double, shape, tol_s0, (double)0.0,
                        mgard_x::error_bound_type::ABS, var_f.data(),
                        compressed_s0, compressed_size_s0, config, false);
    
    mgard_x::decompress(compressed_s0, compressed_size_s0, decompressed_s0, config, false); 

    double rmse_s0   = calc_rmse(var_f.data(), (double *)decompressed_s0, err_s0.data(), num_data); 
    double PE_err_s0 = calc_PE(err_s0.data(), PE_s0.data(), dh, shape[0], shape[1]);
 
    // s=1 compression
    mgard_x::compress(2, mgard_x::data_type::Double, shape, tol_s1, (double)1.0,
                        mgard_x::error_bound_type::ABS, var_f.data(),
                        compressed_s1, compressed_size_s1, config, false);

    mgard_x::decompress(compressed_s1, compressed_size_s1, decompressed_s1, config, false);
    
    double rmse_s1   = calc_rmse(var_f.data(), (double *)decompressed_s1, err_s1.data(), num_data);
    double PE_err_s1 = calc_PE(err_s1.data(), PE_s1.data(), dh, shape[0], shape[1]);

    reader_f.Close();
    
    FILE *fp = fopen((fname_err+"_l2_s0.bin").c_str(), "w");
    fwrite(err_s0.data(), sizeof(double), num_data, fp);
    fclose(fp);
    
    fp = fopen((fname_err+"_l2_s1.bin").c_str(), "w");
    fwrite(err_s1.data(), sizeof(double), num_data, fp);
    fclose(fp);

    fp = fopen((fname_err+"_PE_s0.bin").c_str(), "w");
    fwrite(PE_s0.data(), sizeof(double), num_data, fp);
    fclose(fp);

    fp = fopen((fname_err+"_PE_s1.bin").c_str(), "w");
    fwrite(PE_s1.data(), sizeof(double), num_data, fp);
    fclose(fp);
    
    fp = fopen((fname_err+"_data_ori.bin").c_str(), "w");
    fwrite(var_f.data(), sizeof(double), num_data, fp);
    fclose(fp);

    fp = fopen((fname_err+"_data_s0.bin").c_str(), "w");
    fwrite((double *)decompressed_s0, sizeof(double), num_data, fp);
    fclose(fp);

    fp = fopen((fname_err+"_data_s1.bin").c_str(), "w");
    fwrite((double *)decompressed_s1, sizeof(double), num_data, fp);
    fclose(fp);
    std::cout << "s=0: l2 = " << rmse_s0 << " (rel " << rmse_s0 / (v_max-v_min)<< "), PE_e = " << PE_err_s0 << ", compression ratio = " << (double)(num_data*sizeof(double)) / (double)compressed_size_s0 << "\n";
    std::cout << "s=1: l2 = " << rmse_s1 << " (rel " << rmse_s1 / (v_max-v_min)<< "), PE_e = " << PE_err_s1 << ", compression ratio = " << (double)(num_data*sizeof(double)) / (double)compressed_size_s1 << "\n"; 

    return 0;
}
