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
double potential_energy(T *var_in, T *var_out, std::vector<size_t>dShape)
{
    size_t R=1, C=1, H=1, d=0;
    size_t ndim = dShape.size();
    if (ndim==3) {
        H = dShape[d];
        d++;
    }
    R = dShape[d];
    C = dShape[d+1];

    size_t data_size = R * C * H;
    size_t dim2      = R*C;

//    double average_val = 0.0;
    std::vector<double> var_diff(data_size);
    for (size_t i=0; i<data_size; i++) {
        var_diff[i]  = (double)var_in[i] - (double)var_out[i];
        //average_val += var_diff[i];
    }

    double  PE = 0.0;
    double  ux, uy, uz;
    size_t h_curr, r_curr, k;
    if (ndim==3) {
        std::cout << "h, r, c = " << H << ", " << C << ", " << R << "\n";
        for (size_t h=2; h<H; h++) {
            h_curr = h * dim2;
            for (size_t r=2; r<R; r++) {
                r_curr = r * C;
                for (size_t c=2; c<C; c++) {
                    k   = r_curr + h_curr +c;
                    uz  = (var_diff[k] - var_diff[k-2*dim2])/2.0;
                    uy  = (var_diff[k] - var_diff[k-2*C])/2.0;
                    ux  = (var_diff[k] - var_diff[k-2])/2.0;
                    PE += ux*ux + uy*uy + uz*uz;
                }
            }
        }
    } else if (ndim==2) {
        for (size_t r=2; r<R; r++) {
            r_curr = r * C;
            for (size_t c=2; c<C; c++) {
                k   = r_curr + c;
                uy  = (var_diff[k] - var_diff[k-2*C])/2.0;
                ux  = (var_diff[k] - var_diff[k-2])/2.0;
                PE += ux*ux + uy*uy;
            }
        }
    }
    //std::cout << "averaged data value = " << average_val << "\n";
    return PE;
}

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
    //           1, 2, 3, 4, 5, 6--> rain drop, multiple rains, square, Gaussian sinusoid, plane wave, Gaussian pulse 
    // cfd_cond: 1, 2, 3 --> Dirichlet, Mur boundary condition, Periodic
    int init_fun   = std::stoi(argv[cnt_argv++]);
    int obstacle_t = std::stoi(argv[cnt_argv++]); 
    int cfd_cond   = std::stoi(argv[cnt_argv++]); 
    // simulation space 
    double Dx      = std::stof(argv[cnt_argv++]);
    double Dy      = std::stof(argv[cnt_argv++]);
    // simulation spatial resolution
    double dh      = std::stof(argv[cnt_argv++]);

    size_t Nx = (size_t)std::ceil((double)Dx / dh);
    size_t Ny = (size_t)std::ceil((double)Dy / dh);

    // simulation temporal resolution
    double dt      = std::stof(argv[cnt_argv++]);
    // number of simulation frames
    float T        = std::stof(argv[cnt_argv++]);
    // wave speed
    int n_vp       = std::stoi(argv[cnt_argv++]);
    std::vector<double> wave_c;
    if (n_vp > 0) {
        for (int i=0; i<n_vp; i++) {
            wave_c.push_back(std::stof(argv[cnt_argv++]));
        }
    } else {
        wave_c.resize(Nx*Ny); 
        std::string v_fname(argv[cnt_argv++]);
        FILE *fp = fopen(v_fname.c_str(), "r");
        if (fp == NULL) {
            std::cout << "Error openning velocity mask file\n";
            return -1;
        }  
        fread(wave_c.data(), sizeof(double), Nx*Ny, fp);
        fclose(fp);
    }
    double vp_max = *std::max_element(wave_c.begin(), wave_c.end());
    std::cout << "recommended timestep resolution < " << dh/vp_max/std::sqrt(2) << " based on vp_max = " << vp_max << "\n";
    if (dt >= dh/vp_max/std::sqrt(2)) {
        std::cout << "Need smaller timestep resolution...\n";
        exit(-1);
    }
    // distribution of speed of sound: uniform | gaussian
    std::string material_type(argv[cnt_argv++]);

    // dissipation rate
    double gamma   = std::stof(argv[cnt_argv++]); 
    float temp     = T / dt;
    size_t nframes = (size_t)(temp);
    // relative tolerance. tol = 0 means to write out uncompressed data after each timestep
    double tol_1   = std::stof(argv[cnt_argv++]);
    // using tol_1 for compressing u_curr, and tol_2  for u_next
    double tol_2   = std::stof(argv[cnt_argv++]);
    // error type
    std::string eb_type(argv[cnt_argv++]);
    // write out filename
    std::string fname_wt(argv[cnt_argv++]);
    // readin filename
    std::string fname(argv[cnt_argv++]);
    // initial timestep to read from file
    size_t init_ts = std::stoi(argv[cnt_argv++]);
    // output interval
    size_t wt_interval = (size_t)std::stoi(argv[cnt_argv++]);
    // wave source time
    bool src_on   = std::stoi(argv[cnt_argv++]);    
    double src_ts = std::stof(argv[cnt_argv++]);
    // to write out data after cr_ts
    size_t cr_ts  = std::stoi(argv[cnt_argv++]);

    std::cout << "simulating a domain of [" << Dx << "/" << Nx << ", " << Dy << "/" << Ny << "], at a spacing of " << dh << "\n";
    std::cout << "simulating " << nframes << " steps at a resolution of " << dt << "\n";
    std::cout << "initial function: ";
    if (init_fun==0) std::cout << "read from previous simulation steps in the file: " << fname << " at step " << init_ts << "\n";
    else if (init_fun==1) std::cout << "random raindrop\n";
    else if (init_fun==2) std::cout << "multiple rain drop\n";
    else if (init_fun==3) std::cout << "solid square\n";
    else if (init_fun==4) std::cout << "gaussian sinusoid waves\n";
    else if (init_fun==5) std::cout << "plane waves \n";
    else if (init_fun==6) std::cout << "gaussian pulse source\n";
    
    if (obstacle_t==1) {
        std::cout << "emulating multiple disk obstacles\n";
    }   else if (obstacle_t==2) {
        std::cout << "emulating two slits, recommending using plane wave\n";
    }   else {
        std::cout << "no obstacle\n";
    }
 
    std::cout << "simulating boundary condiction: ";
    if (cfd_cond==1) std::cout << "Dirichlet\n";
    else if (cfd_cond==2) std::cout << "Mur\n";
    else if (cfd_cond==3) std::cout << "Periodic\n";
    std::cout << "dissipation rate = " << gamma << "\n";
    std::cout << "wave speed = ";
    for (int i=0; i<n_vp; i++) std::cout << wave_c[i] << ", ";
    std::cout << "\n";
    if (!strcmp(material_type.c_str(), "uniform")) {
        std::cout << "uniform layered materials\n";    
    } else if (!strcmp(material_type.c_str(), "gaussian")) {
        std::cout << "simulating a bumpped second material\n";
    } else if (!strcmp(material_type.c_str(), "sandwich")) {
        std::cout << "simulating a sandwich material space w/ the middle layer different from the sides\n";
    }
    std::cout << "tolerance = {" << tol_1 << ", " << tol_2 << "} for u_n and u_dt \n";

    if (src_on) std::cout << "turn on a  Gaussian derivative source with a relative ts = " << src_ts << "\n"; 
    
    adios2::ADIOS ad(MPI_COMM_WORLD);
    adios2::IO writer_io = ad.DeclareIO("Output");

    adios2::Engine writer = writer_io.Open(fname_wt, adios2::Mode::Write);
    adios2::Variable<double> variable_u, variable_u_dt, variable_cr;
    variable_u  = writer_io.DefineVariable<double>("u_data" , adios2::Dims{Nx, Ny}, adios2::Dims{0,0}, adios2::Dims{Nx, Ny});
    std::vector<double> u_dt, u_at;
    if (tol_1>0 || tol_2>0) {
        // compression ratio of u_n and u_n-u_np1
        variable_cr = writer_io.DefineVariable<double>("u_CR"   , adios2::Dims{2}, adios2::Dims{0}, adios2::Dims{2});
        variable_u_dt = writer_io.DefineVariable<double>("u_dt", adios2::Dims{Nx, Ny}, adios2::Dims{0,0}, adios2::Dims{Nx, Ny});
        u_dt.resize(Nx*Ny);
        u_at.resize(Nx*Ny);
    }
    
    double max_intensity = 10.0;
    // for Gaussian wave
    double freq = 2.0 * M_PI / T;
    // for rain drops
    float drop_probability = 1;
    std::vector<double> gauss_template;
    double NDx = 12, NDy=24;
    size_t n_drops = 4;
    // for Gaussian pulse source
    double f0 = 100; // dominant frequency of the source (Hz)
    double t0 = 0.1; // source time shift (s) 
    size_t srcx = size_t(Dx * 0.5 / dh); 
    size_t srcy = size_t(Dy * 0.5 / dh);
    size_t src_pos = srcx * Ny + srcy;
    double src_intensity = 1.0;
    // for compression
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
    WaveEquation <double> waveSim(Nx-1, Ny-1, 0, dt, dh, gamma, cfd_cond); 
    
    // initialize wave velocity
    std::vector<double> speed_sound(Nx*Ny);
    if (!strcmp(material_type.c_str(), "uniform")) {
        velocity_Layered_uniform(speed_sound.data(), wave_c, n_vp, Nx, Ny, 1);
    } else if (!strcmp(material_type.c_str(), "gaussian")) {
        std::fill(speed_sound.begin(), speed_sound.end(), wave_c[0]);
        double peak_h = 0.5;
        velocity_Gaussian_2d(speed_sound.data(), wave_c[1], peak_h, Nx, Ny);
    } else if (!strcmp(material_type.c_str(), "sandwich")) {
        std::fill(speed_sound.begin(), speed_sound.end(), wave_c[0]);
        size_t width_slit = (size_t)(0.04 * (double)Ny);
        bool flg_v = velocity_sandwich(speed_sound.data(), wave_c[1], width_slit, Nx, Ny, 1);
        if (flg_v==false) {
            std::cout << "Fail on creating a sandwich material space\n";
            exit(-1);
        }
    } else if (!strcmp(material_type.c_str(), "mask")) {
        velocity_mask(speed_sound.data(), wave_c.data(), Nx, Ny, 1);
    }      
    
    //FILE *fp = fopen("velocity_sandwich.bin", "w");
    //fwrite(speed_sound.data(), sizeof(double), Nx*Ny, fp);
    //fclose(fp);
    std::cout << "velocity = " << *std::max_element(speed_sound.begin(), speed_sound.end()) << "\n";
    waveSim.init_vp(speed_sound.data());
    double *obstacle_m = NULL;
    if (obstacle_t) {
        switch (obstacle_t) {
            case 1: {
                size_t n_disks    = 8;
                size_t max_radius = size_t(0.015 * Nx);  
                std::vector<size_t> x_pos = {311, 383, 221, 19, 272, 235, 571, 477, 430, 457};
                std::vector<size_t> y_pos = {426, 182, 223, 37, 9, 94, 148, 29, 245, 295};
                obstacle_m        = emulate_N_disk<double>(Nx, Ny, n_disks, max_radius, x_pos, y_pos);
                break;
            }
            case 2: {
                size_t width        = size_t ((1.0/16.0) * Ny);
                size_t p1           = size_t ((5.0/16.0) * Ny);
                size_t p2           = size_t ((5.0/8.0 ) * Ny);
                size_t thickness    = size_t ((1.0/32.0) * Nx);
                size_t dist_to_edge = size_t ((1.0/5.0) * Nx);
                obstacle_m = emulate_two_slits<double>(Nx, Ny, p1, p2, width, dist_to_edge, thickness);
                break;
            }
        }
    }

    if ((init_fun==1) || (init_fun==2)) {
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
                if ((r-cx)*(r-cx)+(c-cy)*(c-cy) < NDx*NDx) {
                    gauss_template[r*NDy+c] = max_intensity * exp(-(px[r]*px[r] + py[c]*py[c]));
                }
            }
        }
    }
    switch (init_fun) {
        case 0: {
            adios2::IO reader_io = ad.DeclareIO("Input");
            reader_io.SetEngine("BP");
            adios2::Engine reader = reader_io.Open(fname, adios2::Mode::ReadRandomAccess);
            adios2::Variable<double> variable;
            if ((tol_1==0) && (tol_2==0) && (init_ts>0)) { // optimized checkpoint restart 
                // load the checkpoint data compressed using optimized approach 
                std::cout << "load u_dt\n";
                variable = reader_io.InquireVariable<double>("u_dt"); 
            } else {
                // load the checkpoint data compressed using non-optimized approach
                std::cout << "load u_{n}\n";
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
            std::cout << "Previous step " << init_ts << ", now at step " << reader.CurrentStep() << "\n";

            std::cout << "begin to load the current step...\n";
            if ((init_ts>0)  && (tol_1==0) && (tol_2==0)) {  
                std::cout << "read u_at\n";
                variable = reader_io.InquireVariable<double>("u_data");
                variable.SetStepSelection({init_ts, 1});
            } else {
                std::cout << "read u_{n+1}\n";
                variable.SetStepSelection({init_ts+1, 1});
            }
            reader.Get(variable, init_v.data());
            reader.PerformGets();
            std::cout <<  "u_curr: min/max = " << *std::min_element(init_v.begin(), init_v.end()) << ", " << *std::max_element(init_v.begin(), init_v.end()) << "\n";
            if ((init_ts>0) && (tol_1==0) && (tol_2==0)) {  // change (u_n - u_nm1)/2 and (u_n+u_nm1)/2 to u_n and u_mn1
                std::transform(init_vpre.begin(), init_vpre.end(), init_vpre.begin(), [](double num) { return num / 2.0; });
                std::vector<double> u_A(init_v);
                // u_n = u_A + 0.5*u_dt
                std::transform(init_v.begin(), init_v.end(), init_vpre.begin(), init_v.begin(), std::plus<double>()); 
                // u_{n-1} = u_A - 0.5*u_dt
                std::transform(u_A.begin(), u_A.end(), init_vpre.begin(), init_vpre.begin(), std::minus<double>());
                u_A.clear();
            }
            waveSim.init_u_n(init_vpre.data());
            waveSim.init_u_np1(init_v.data());
            reader.Close();
            init_v.clear();
            init_vpre.clear(); 
            if ((tol_1==0) && (tol_2==0)) {
                writer.BeginStep();
                writer.Put<double>(variable_u, waveSim.u_n.data(), adios2::Mode::Sync);
                writer.PerformPuts();
                writer.EndStep();
                //writer.BeginStep();
                //writer.Put<double>(variable_u, waveSim.u_np1.data(), adios2::Mode::Sync);
                //writer.PerformPuts();
                //writer.EndStep();
            }
            break;
        }
        case 1: { 
            fun_rainDrop<double>(waveSim.u_np1.data(), Nx, Ny, 1, NDx, NDy, 1, gauss_template.data(), drop_probability);
            break;
        }
        case 2: {
            fun_MultiRainDrop<double>(waveSim.u_np1.data(), Nx, Ny, 1, NDx, NDy, 1, gauss_template.data(), drop_probability, n_drops);
            break;
        }
        case 3: {
            fun_square<double>(waveSim.u_np1.data(), Nx, Ny, 1, NDx, NDy, 1, max_intensity);
            break;
        }
        case 4: {
            int n_waves = 4;
            std::vector<double> sigma = {50, 25, 35, 60};
            std::vector<double> intensity = {25, 30, 45, 15};
            std::vector<double> freq_x = {freq*T/18.0, freq*T/18.0, freq*T/24.0, freq*T/22.0};
            std::vector<double> freq_y = {freq*T/18.0, freq*T/32.0, freq*T/17.0, freq*T/18.0};
            std::vector<double> freq_z(n_waves, 0.0);
            fun_gaussian_wave(waveSim.u_np1.data(), Nx, Ny, 1, intensity, sigma, freq_x, freq_y, freq_z, n_waves);
            break;
        }
        case 5: {
            size_t n_waves = size_t(0.02 * (double)Nx);
            size_t pos_x   = size_t(0.3  * (double)Nx); 
            fun_plane_waves<double>(waveSim.u_np1.data(), Nx, Ny, pos_x, 5*max_intensity, freq*T/20.0, n_waves);
            break;
        }
        //case 6: {
        //    fun_Gaussian_pulse(waveSim.u_np1.data(), f0, t0, src_intensity, srcx, srcy, (size_t)0, Ny, (size_t)1); 
        //}
        default:
            break; 
    }
    if (obstacle_t) { 
            // compression & decompression will introduce "noise", changing the obstacle occupied space to non-zero
            // apply Dirichlet boundary condition to the scattering on obstacles
            std::transform(waveSim.u_np1.begin(), waveSim.u_np1.end(), obstacle_m, waveSim.u_np1.begin(), std::multiplies<double>());
            std::transform(waveSim.u_n.begin(), waveSim.u_n.end(), obstacle_m, waveSim.u_n.begin(), std::multiplies<double>());
    }

    drop_probability  = 0.01; 
    size_t iter_frame = 0;
    while (iter_frame<nframes) {
        if (iter_frame % wt_interval == 0) std::cout << iter_frame << "/" << nframes << "\n";
        // source update
        if ((src_on) && (init_fun || iter_frame || ((tol_1*tol_2==0) && (init_ts)>0))) {
            double ts = (init_fun==0) ? ((tol_1*tol_2==0) ? 2 : 1) : 0; /* c.r. from t=2*/;
            ts += (double)(iter_frame + 1) + src_ts /* secondary c.r. */ + init_ts;
            t0 = (ts>=1875) ? dt*2000 : 0.1;
            if ((ts<125) || ((ts>=1875) && (ts<2125))) {
                waveSim.u_np1[src_pos] = src_Gaussian_pulse(f0, ts*dt-t0, src_intensity);
                printf("update pos %ld, src[%ld]=%.5e\n",src_pos, (size_t)ts, waveSim.u_np1[src_pos]); 
            }
        }
        // wave update
        waveSim.update_2d();

        if (obstacle_t) {
            // apply Dirichlet boundary condition to the scattering on obstacles
            std::transform(waveSim.u_np1.begin(), waveSim.u_np1.end(), obstacle_m, waveSim.u_np1.begin(), std::multiplies<double>());
        }
        if ((iter_frame % wt_interval == 0) && (iter_frame>=cr_ts)) {
            writer.BeginStep();
            if (tol_1>0 && tol_2>0) {
                void *compressed_array_cpu  = NULL;
                void *compressed_array2_cpu = NULL;
                std::transform(waveSim.u_np1.begin(), waveSim.u_np1.end(), waveSim.u_n.begin(), u_at.begin(), std::plus<double>());
                std::transform(u_at.begin(), u_at.end(), u_at.begin(), [](double num) { return num / 2.0; });
                std::transform(waveSim.u_np1.begin(), waveSim.u_np1.end(), waveSim.u_n.begin(), u_dt.begin(), std::minus<double>());
                if (obstacle_t) {
                    // apply Dirichlet boundary condition to the scattering on obstacles
                    std::transform(u_dt.begin(), u_dt.end(), obstacle_m, u_dt.begin(), std::multiplies<double>());
                    std::transform(u_at.begin(), u_at.end(), obstacle_m, u_at.begin(), std::multiplies<double>());
                }
                if (strcmp(eb_type.c_str(), "ABS")==0) {
                    mgard_x::compress(2, mgard_x::data_type::Double, shape, tol_1, s1,
                        mgard_x::error_bound_type::ABS, u_at.data(),
                        compressed_array_cpu, compressed_size_u, config, false);
                
                    mgard_x::compress(2, mgard_x::data_type::Double, shape, tol_2, s2,
                        mgard_x::error_bound_type::ABS, u_dt.data(),
                        compressed_array2_cpu, compressed_size_dt, config, false);
                } else {
                    mgard_x::compress(2, mgard_x::data_type::Double, shape, tol_1, s1,
                        mgard_x::error_bound_type::REL, u_at.data(),
                        compressed_array_cpu, compressed_size_u, config, false);

                    mgard_x::compress(2, mgard_x::data_type::Double, shape, tol_2, s2,
                        mgard_x::error_bound_type::REL, u_dt.data(),
                        compressed_array2_cpu, compressed_size_dt, config, false);
                }
                compression_ratio[0] = data_bytes / (double) compressed_size_u;
                compression_ratio[1] = data_bytes / (double) compressed_size_dt;
                std::cout << "compression ratios = " << compression_ratio[0] << " and " << compression_ratio[1] << ", averaged = " << data_bytes*2 / (double) (compressed_size_u + compressed_size_dt) << "\n";
                void *decompressed_array_cpu = NULL;
                void *decompressed_array2_cpu = NULL;
                mgard_x::decompress(compressed_array_cpu, compressed_size_u,
                    decompressed_array_cpu, config, false);
                mgard_x::decompress(compressed_array2_cpu, compressed_size_dt,
                    decompressed_array2_cpu, config, false);
            
                //error_calc(waveSim.u_np1, (double *)decompressed_array_cpu, tol_1); 
                //error_calc(u_dt, (double *)decompressed_array2_cpu, tol_2);
                writer.Put<double>(variable_u   , (double *)decompressed_array_cpu , adios2::Mode::Sync);
                writer.Put<double>(variable_u_dt, (double *)decompressed_array2_cpu, adios2::Mode::Sync);
                writer.Put<double>(variable_cr, compression_ratio.data(), adios2::Mode::Sync);
            } else if (tol_1>0 && tol_2==0) {  // non-optimized compression
                void *compressed_array_cpu  = NULL;
                if (strcmp(eb_type.c_str(), "ABS")==0) {
                    mgard_x::compress(2, mgard_x::data_type::Double, shape, tol_1, s2,
                        mgard_x::error_bound_type::ABS, waveSim.u_n.data(),
                        compressed_array_cpu, compressed_size_u, config, false);
                } else {
                    mgard_x::compress(2, mgard_x::data_type::Double, shape, tol_1, s2,
                        mgard_x::error_bound_type::REL, waveSim.u_n.data(),
                        compressed_array_cpu, compressed_size_u, config, false);
                }
                std::cout << "compression ratios = " << data_bytes / (double) compressed_size_u << "\n";
                void *decompressed_array_cpu = NULL;
                mgard_x::decompress(compressed_array_cpu, compressed_size_u,
                    decompressed_array_cpu, config, false);
                writer.Put<double>(variable_u   , (double *)decompressed_array_cpu , adios2::Mode::Sync);
                //std::vector<size_t> dShape = {Nx, Ny};
                //double PE = potential_energy(waveSim.u_n.data(), (double *)decompressed_array_cpu, dShape);
                //printf("PE = %.8f, eb / sqrt(PE) = %.8f\n", std::sqrt(PE), tol_1 / std::sqrt(PE));
            } 
            else { // no compression
                writer.Put<double>(variable_u, waveSim.u_n.data(), adios2::Mode::Sync);
            }
            writer.PerformPuts();
            writer.EndStep();
        }
        iter_frame ++;
    }
    writer.Close();
    
    u_dt.clear();
    u_at.clear();
    delete [] obstacle_m;
    MPI_Finalize();
    return 0;
}













