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

double calc_KE(double *u_prev, double *u, double *velocity, double dt, size_t num_data)
{
    double diff_t, KE = 0.0;
    double dt2 = dt*dt;
    for (size_t i=0; i<num_data; i++) {
        diff_t = (u[i] - u_prev[i]);
        KE += diff_t * diff_t / dt2 / velocity[i];
    }
    return KE * 0.5;
}

double calc_PE(double *u, double dh, size_t Nx, size_t Ny)
{
    double PE = 0.0;
    double ux, uy;
    size_t r_curr, k;
    for (size_t r=1; r<Nx; r++) {
        r_curr = r * Ny;
        for (size_t c=1; c<Ny; c++) {
            k   = r_curr + c;
            ux  = (u[k] - u[k-Ny])/dh;
            uy  = (u[k] - u[k-1])/dh;
            PE += ux*ux + uy*uy;
        }
    }
    return PE * 0.5; 
}

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
    std::string fname_g(argv[cnt_argv++]);
    // velocity data
    std::string fname_v(argv[cnt_argv++]);
    double C;
    if (strcmp(fname_v.c_str(), "NaN")==0) {
        // constant velocity across the space
        C = std::stof(argv[cnt_argv++]);
    }
    // simulation spatial resolution
    double dh      = std::stof(argv[cnt_argv++]);
    // simulation temporal resolution
    double dt      = std::stof(argv[cnt_argv++]);
    // from the init_ts step in frame_f to compare against frame_g 
    size_t init_ts = std::stoi(argv[cnt_argv++]);    
    std::string fname_err(argv[cnt_argv++]);
    size_t total_Steps = std::stoi(argv[cnt_argv++]);

    std::cout << "original file: " << fname_f.c_str() << "\n";
    std::cout << "c.r. file: " << fname_g.c_str() << "\n";
    std::cout << "velocity mask file: " << fname_v.c_str() << "\n";
    std::cout << "output file: " << fname_err.c_str() << "\n"; 
    std::cout << "dt = " << dt << ", dh = " << dh << ", init_ts = " << init_ts << "\n";

    adios2::ADIOS ad;
    adios2::IO reader_io_f = ad.DeclareIO("Original");
    adios2::IO reader_io_g = ad.DeclareIO("Lossy");
    reader_io_f.SetEngine("BP");
    reader_io_g.SetEngine("BP");
    adios2::Engine reader_f = reader_io_f.Open(fname_f, adios2::Mode::ReadRandomAccess);
    adios2::Engine reader_g = reader_io_g.Open(fname_g, adios2::Mode::ReadRandomAccess);

    adios2::Variable<double> variable_f = reader_io_f.InquireVariable<double>("u_data");
    adios2::Variable<double> variable_g = reader_io_g.InquireVariable<double>("u_data");
    size_t available_Steps = std::min(variable_g.Steps(), variable_f.Steps() - init_ts); 
    if (total_Steps>0) {
        total_Steps = std::min(total_Steps, available_Steps); 
    }
    std::cout << "total number of steps: " << variable_f.Steps() << ", read from " << init_ts << " to " << total_Steps << " timestep \n";
    reader_g.Close();
    reader_g = reader_io_g.Open(fname_g, adios2::Mode::Read);
 
    std::vector<std::size_t> shape = variable_f.Shape();
    size_t num_data = shape[0]*shape[1];
    std::vector<double> var_f(num_data);
    std::vector<double> var_g(num_data);
    // difference data
    std::vector<double> var_e(num_data);
    std::vector<double> var_prev(num_data);
    std::vector<double> rmse(total_Steps);
    std::vector<double> PE_e(total_Steps);
    std::vector<double> KE_e(total_Steps);

    std::vector<double> wave_c(num_data);
    if (strcmp(fname_v.c_str(), "NaN")==0) {
        std::fill(wave_c.begin(), wave_c.end(), C); 
    } else {
        FILE *fp = fopen(fname_v.c_str(), "r");
        if (fp == NULL) {
            std::cout << "Error openning velocity mask file\n";
            return -1;
        }
        fread(wave_c.data(), sizeof(double), num_data, fp);
        fclose(fp);
    }

    size_t cnt = 0;
    while (true) {
        // Begin step
        adios2::StepStatus read_status = reader_g.BeginStep(adios2::StepMode::Read, 10.0f);
        if (read_status == adios2::StepStatus::NotReady) {
            // std::cout << "Stream not ready yet. Waiting...\n";
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            continue;
        }
        else if (read_status != adios2::StepStatus::OK) {
            break;
        }
        variable_f.SetStepSelection({init_ts+cnt, 1});
        reader_f.Get(variable_f, var_f);
        reader_g.Get(variable_g, var_g); 
        reader_f.PerformGets();
        reader_g.PerformGets();
         
        rmse[cnt] = calc_rmse(var_f.data(), var_g.data(), var_e.data(), num_data);
        if (cnt>0) {
            KE_e[cnt] = calc_KE(var_prev.data(), var_e.data(), wave_c.data(), dt, num_data);
        }
        PE_e[cnt] = calc_PE(var_e.data(), dh, shape[0], shape[1]);
        std::copy(var_e.begin(), var_e.end(), var_prev.begin());
        reader_g.EndStep(); 
        if (cnt % 500 == 0) std::cout << cnt << " / " << init_ts+cnt << "\n"; 
        cnt ++;
        if (cnt == total_Steps) break;
    }
    reader_f.Close();
    reader_g.Close();
    
    FILE *fp = fopen((fname_err+"_l2.bin").c_str(), "w");
    fwrite(rmse.data(), sizeof(double), num_data, fp);
    fclose(fp);
    
    fp = fopen((fname_err+"_KE_e.bin").c_str(), "w");
    fwrite(KE_e.data(), sizeof(double), num_data, fp);
    fclose(fp);

    fp = fopen((fname_err+"_PE_e.bin").c_str(), "w");
    fwrite(PE_e.data(), sizeof(double), num_data, fp);
    fclose(fp);

    for (size_t i=0; i<500; i++) {
        std::cout << i << ": l2 = " << rmse[i] << ", KE_e = " << KE_e[i] << ", PE_e = " << PE_e[i] << "\n";
    }
    return 0;
}
