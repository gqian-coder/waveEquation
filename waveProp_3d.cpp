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

// center_x, center_y: percentage, location in the x and y-axis
// A: magnitude of the sinusoid signal
// iter: iteration of the update
// freq: frequency -- 2*M_PI/T
template <typename Real>
void fun_sinusoidal(Real *u, size_t Nx, size_t Ny, size_t Nz, float center_x, float center_y, 
                    float center_z, Real A, size_t iter, Real freq)
{
    if ((center_x <= 0) || (center_x>=1) || (center_y<=0) || (center_y>=1) || (center_z<=0) || (center_z>=1)) {
        std::cout << "The initial position of the sinusoidal ring must be a positive number and less than the boundary\n";
        return;
    }
    size_t pos_x = static_cast<size_t>(Nx*center_x);
    size_t pos_y = static_cast<size_t>(Ny*center_y);
    size_t pos_z = static_cast<size_t>(Nz*center_z);
    u[pos_x * Ny * Nz + pos_y * Nz + pos_z] = A * sin(iter*freq); 
}

// A: magnitude
template <typename Real>
void fun_velocity(Real *u, size_t Nx, size_t Ny, size_t Nz, 
                  Real dx, Real dy, Real dz, Real A, Real sigma)
{
    size_t offset_x, offset_y, k;
    std::vector <Real> px(Nx), py(Ny), pz(Nz); 
    for (size_t r=0; r<Nx; r++) {
        px[r] = (static_cast<Real>(r)-Nx)*dx;
    }
    for (size_t c=0; c<Ny; c++) {
        py[c] = (static_cast<Real>(c)-Ny)*dy;
    }
    for (size_t z=0; z<Nz; z++) {
        pz[z] = (static_cast<Real>(z)-Nz)*dz;
    }
    for (size_t r=0; r<Nx; r++) {
        offset_x = r * Ny * Nz;
        for (size_t c=0; c<Ny; c++) {
            offset_y = offset_x + c * Nz;
            for (size_t z=0; z<Nz; z++) {
                k  = offset_y + z;
                u[k] = A * exp(-((px[r]*px[r] + py[c]*py[c] + pz[z]*pz[z])/sigma)); 
            }
        }
    } 
}

// initialization of a solid square
// NDx, NDy: square's dimension
template<typename Real>
void fun_square(Real *u, size_t Nx, size_t Ny, size_t Nz, 
                size_t NDx, size_t NDy, size_t NDz, Real intensity)
{
    std::cout << "fun_square\n";
    size_t cx = (size_t)(NDx/2);
    size_t cy = (size_t)(NDy/2);
    size_t cz = (size_t)(NDz/2);
    std::random_device rd;  // Non-deterministic random device (for seeding)
    std::mt19937 gen(rd()); // Mersenne Twister pseudo-random generator
    std::uniform_real_distribution<> dis(0.0, 1.0);
    float random_x = dis(gen);
    float random_y = dis(gen);
    float random_z = dis(gen);
    size_t x = static_cast<size_t>(random_x * (Nx-cx-1) + cx);
    size_t y = static_cast<size_t>(random_y * (Ny-cy-1) + cy);
    size_t z = static_cast<size_t>(random_z * (Nz-cz-1) + cz);
    for (size_t r=0; r<NDx; r++) {
        for (size_t c=0; c<NDy; c++) {
            for (size_t h=0; h<NDz; h++) {
                u[(x-cx + r)*Ny*Nz + (y-cy+c)*Nz + h-cz+z] = intensity; 
                //std::cout << gauss_template[r*NDy+c] << "\n";
            }
        }
    }
}


// drop_probability: Raindrop probability (with each time tick) and intensity.
// NDx, NDy: droplet's region
// one droplet at each time
template <typename Real>
void fun_rainDrop(Real *u, size_t Nx, size_t Ny, size_t Nz, 
                  size_t NDx, size_t NDy, size_t NDz, 
                  Real *gauss_template, float drop_probability)
{
    size_t cx = (size_t)(NDx/2);
    size_t cy = (size_t)(NDy/2);
    size_t cz = (size_t)(NDz/2);
    std::random_device rd;  // Non-deterministic random device (for seeding)
    std::mt19937 gen(rd()); // Mersenne Twister pseudo-random generator
    std::uniform_real_distribution<> dis(0.0, 1.0);
    float random_number = dis(gen); 
    //std::cout << "probability: " << random_number << "\n"; 
    if (random_number < drop_probability) {
        float random_x = 0.5;//dis(gen);
        float random_y = 0.5;//dis(gen);
        float random_z = 0.5;//dis(gen);
        size_t x = static_cast<size_t>(random_x * (Nx-cx-1) + cx);
        size_t y = static_cast<size_t>(random_y * (Ny-cy-1) + cy);
        size_t z = static_cast<size_t>(random_z * (Nz-cz-1) + cz);
        for (size_t r=0; r<NDx; r++) {
            for (size_t c=0; c<NDy; c++) { 
                for (size_t h=0; h<NDz; h++) {
                    u[(x-cx + r)*Ny*Nz + (y-cy+c)*Nz + h-cz+z] = gauss_template[r*NDy*NDz+c*NDz+h]; 
                }
            }            
        }
    }
    //std::cout << *std::min_element(u, u+Nx*Ny*Nz) << ", " << *std::max_element(u, u+Nx*Ny*Nz) << "\n";
}

// drop_probability: Raindrop probability (with each time tick) and intensity.
// drop multiple droplets
template <typename Real>
void fun_MultiRainDrop(Real *u, size_t Nx, size_t Ny, size_t Nz, 
                       size_t NDx, size_t NDy, size_t NDz,
                       Real *gauss_template, float drop_probability, int nDrops)
{
    size_t cx = (size_t)(NDx/2);
    size_t cy = (size_t)(NDy/2);
    size_t cz = (size_t)(NDz/2);
    std::random_device rd;  // Non-deterministic random device (for seeding)
    std::mt19937 gen(rd()); // Mersenne Twister pseudo-random generator
    std::uniform_real_distribution<> dis(0.0, 1.0);
    float random_number = dis(gen);
    //std::cout << "probability: " << random_number << "\n";
    for (int d=0; d<nDrops; d++) {     
        if (random_number < drop_probability) {
            float random_x = dis(gen);
            float random_y = dis(gen);
            float random_z = dis(gen);
            size_t x = static_cast<size_t>(random_x * (Nx-cx-1) + cx);
            size_t y = static_cast<size_t>(random_y * (Ny-cy-1) + cy);
            size_t z = static_cast<size_t>(random_z * (Nz-cz-1) + cz);
            std::cout << "x, y, z = " << x << ", " << y  << ", " << z << "\n";
            for (size_t r=0; r<NDx; r++) {
                for (size_t c=0; c<NDy; c++) {
                    for (size_t h=0; h<NDz; h++) {
                        float intensity = dis(gen);
                        u[(x-cx + r)*Ny*Nz + (y-cy+c)*Nz + h-cz+z] = intensity*gauss_template[r*NDy*NDz+c*NDz+h];
                    }
                }
            }
        }
    }
    //std::cout << *std::min_element(u, u+Nx*Ny) << ", " << *std::max_element(u, u+Nx*Ny) << "\n";
}



int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, np_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np_size);

    int cnt_argv = 1;
   
    // init_fun: 0 --> read from a data file
    //           1, 2, 3, 4, 5 --> sinusoidal, rain drop, multiple rains, square, velocity 
    // cfd_cond: 1, 2 --> Dirichlet, Mur boundary condition
    int init_fun   = std::stoi(argv[cnt_argv++]);
    bool continuous_update; 
    std::istringstream(argv[cnt_argv++]) >> std::boolalpha >> continuous_update;
    bool iterative_sim;
    std::istringstream(argv[cnt_argv++]) >> std::boolalpha >> iterative_sim;
    int cfd_cond   = std::stoi(argv[cnt_argv++]); 
    // simulation space 
    double Dx      = std::stof(argv[cnt_argv++]);
    double Dy      = std::stof(argv[cnt_argv++]);
    double Dz      = std::stof(argv[cnt_argv++]);
    // simulation spatial resolution
    double dh      = std::stof(argv[cnt_argv++]);
    // simulation temporal resolution
    double dt      = std::stof(argv[cnt_argv++]);
    // number of simulation frames
    float T        = std::stof(argv[cnt_argv++]);
    float temp     = T / dt;
    // wave speed
    double C       = std::stof(argv[cnt_argv++]);
    // dissipation rate
    double gamma   = std::stof(argv[cnt_argv++]);    
    size_t nframes = (size_t)(temp);
    // relative tolerance. tol = 0 means to write out uncompressed data after each timestep
    double tol     = std::stof(argv[cnt_argv++]);
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
    size_t Nz = (size_t)std::ceil((double)Dz / dh);

    std::cout << "simulating a domain of [" << Dx << "/" << Nx << ", " << Dy << "/" << Ny << ", " << Dz << "/" << Nz << "], at a spacing of " << dh << "\n"; 
    std::cout << "simulating " << nframes << " steps at a resolution of " << dt << "\n";
    std::cout << "initial function: ";
    if (init_fun==0) std::cout << "read from previous simulation steps in the file: " << fname << "\n";
    else if (init_fun==1) std::cout << "sinusoidal \n";
    else if (init_fun==2) std::cout << "random raindrop\n";
    else if (init_fun==3) std::cout << "multiple rain drop\n";
    else if (init_fun==4) std::cout << "solid square\n";
    else std::cout << "exponential wave velocity\n";
    std::cout << "simulating boundary condiction: ";
    if (cfd_cond==1) std::cout << "Dirichlet\n";
    else std::cout << "Mur\n"; 

    adios2::ADIOS ad(MPI_COMM_WORLD);
    adios2::IO writer_io = ad.DeclareIO("Output");

    adios2::Engine writer = writer_io.Open(fname_wt, adios2::Mode::Write);
    adios2::Variable<double> variable_wt;
    variable_wt = writer_io.DefineVariable<double>("u_data", adios2::Dims{Nx, Ny, Nz}, adios2::Dims{0,0,0}, adios2::Dims{Nx, Ny, Nz});
    /*
    adios2::Operator op = ad.DefineOperator("mgard", "mgard");
    if (tol>0){
        std::cout << "save the compressed data per timestep w/ eb = " << tol << "\n";
        variable_wt.AddOperation(op, {{"accuracy", std::to_string(tol)}, {"mode", "ABS"}}); 
    }
    */
    float max_intensity = 10.0;
    float center_x = 0.5, center_y = 0.5, center_z = 0.5;
    float freq = 2.0 * M_PI / T;
    float drop_probability = 1;
    std::vector<double> gauss_template;
    double NDx = 12, NDy=24, NDz=12;
    size_t n_drops = 5;

    // compression parameters
    mgard_x::Config config;
    config.lossless = mgard_x::lossless_type::Huffman_Zstd;
    config.dev_type = mgard_x::device_type::CUDA;
    //config.dev_id   = 1;
    std::vector<mgard_x::SIZE> shape{Nx, Ny, Nz};
    double s = 0.0;
    size_t compressed_size_step = 0;
    WaveEquation <double> waveSim(Nx-1, Ny-1, Nz-1, dt, dh, C, gamma, cfd_cond); 

    if ((init_fun==2) || (init_fun==3)) {
        // Width of the Gaussian profile for each initial drop.
        double drop_width = 5e-2 * Dx;
        // Size of the Gaussian template each drop is based on.
        NDx = NDy = NDz = (size_t) std::ceil(drop_width / dh);
        size_t cx = (size_t)(NDx/2);
        size_t cy = (size_t)(NDy/2);
        size_t cz = (size_t)(NDz/2);
        gauss_template.resize(NDx*NDy*NDz);
        std::vector <double> px(NDx), py(NDy), pz(NDz);
        for (size_t r=0; r<NDx; r++) {
            px[r] = ((double)r - cx)/drop_width;
        }
        for (size_t c=0; c<NDy; c++) {
            py[c] = (double(c)-cy)/drop_width;
        }
        for (size_t z=0; z<NDz; z++) {
            pz[z] = (double(z)-cz)/drop_width;
        }
        std::cout << "rainDrop region: " << NDx << " x " << NDy << " x " << NDz <<" pixels\n";
        for (size_t r=0; r<NDx; r++) {
            size_t offset_r = r * NDy * NDz;
            for (size_t c=0; c<NDy; c++) {
                size_t offset_c = offset_r + c * NDz;
                for (size_t z=0; z<NDz; z++) {
                    gauss_template[offset_c+z] = max_intensity * exp(-(px[r]*px[r] + py[c]*py[c] + pz[z]*pz[z]));
                }
            }
        }
    }
    switch (init_fun) {
        case 0: {
            adios2::IO reader_io = ad.DeclareIO("Input");
            reader_io.SetEngine("BP");
            adios2::Engine reader = reader_io.Open(fname, adios2::Mode::ReadRandomAccess);
            adios2::Variable<double> variable = reader_io.InquireVariable<double>("u_data");
            std::cout << "total number of steps: " << variable.Steps() << ", read from " << init_ts << " timestep \n";
            double* init_v = (double *)malloc(sizeof(double) * Nx*Ny*Nz);
            variable.SetSelection({adios2::Dims{0,0,0}, adios2::Dims{Nx, Ny, Nz}});
            variable.SetStepSelection({init_ts, 1}); 
            reader.Get(variable, init_v);
            reader.PerformGets();
            std::cout << *std::min_element(init_v, init_v+Nx*Ny*Nz) << ", " << *std::max_element(init_v, init_v+Nx*Ny*Nz) << "\n";
            waveSim.init_u_n(init_v); 
            variable.SetStepSelection({init_ts+1, 1});
            reader.Get(variable, init_v);
            reader.PerformGets();
            std::cout << *std::min_element(init_v, init_v+Nx*Ny*Nz) << ", " << *std::max_element(init_v, init_v+Nx*Ny*Nz) << "\n";
            waveSim.init_u_np1(init_v);
            reader.Close();
            break;
        }
        case 1:
            break;
        case 2: { 
            fun_rainDrop<double>(waveSim.u_np1.data(), Nx, Ny, Nz, NDx, NDy, NDz, gauss_template.data(), drop_probability);
            break;
        }
        case 3: {
            fun_MultiRainDrop<double>(waveSim.u_np1.data(), Nx, Ny, Nz, NDx, NDy, NDz, gauss_template.data(), drop_probability, n_drops);
            break;
        }
        case 4: {
            fun_square<double>(waveSim.u_np1.data(), Nx, Ny, Nz, NDx, NDy, NDz, max_intensity);
            break;
        }
        case 5: {
            double sigma = 0.1;
            fun_velocity<double>(waveSim.u_np1.data(), Nx, Ny, Nz, dh, dh, dh, max_intensity, sigma);
            break;
        }
        default:
            break; 
    }

    drop_probability  = 0.01; 
    size_t iter_frame = 0;
    while (iter_frame<nframes) {
        std::cout << iter_frame << "/" << nframes << "\n";
        writer.BeginStep();
        if (continuous_update) {
            switch (init_fun) {
                case 1:
                    fun_sinusoidal<double>(waveSim.u_np1.data(), Nx, Ny, Nz, center_x, center_y, center_z, max_intensity, iter_frame, freq);
                    break;
                case 2:
                    fun_rainDrop<double>(waveSim.u_np1.data(), Nx, Ny, Nz, NDx, NDy, NDz, gauss_template.data(), drop_probability);
                    break;
                case 3:
                    fun_MultiRainDrop<double>(waveSim.u_np1.data(), Nx, Ny, Nz, NDx, NDy, NDz, gauss_template.data(), drop_probability, n_drops);
                default:
                    break;
            }
        }
        waveSim.update_3d(); 
        if (tol>0) {
            void *compressed_array_cpu = NULL;
            if (strcmp(eb_type.c_str(), "ABS")==0) {
                mgard_x::compress(3, mgard_x::data_type::Double, shape, tol, s,
                    mgard_x::error_bound_type::ABS, waveSim.u_np1.data(),
                    compressed_array_cpu, compressed_size_step, config, false);
            } else {
                mgard_x::compress(3, mgard_x::data_type::Double, shape, tol, s,
                    mgard_x::error_bound_type::REL, waveSim.u_np1.data(),
                    compressed_array_cpu, compressed_size_step, config, false);
            }
            void *decompressed_array_cpu = NULL;
            mgard_x::decompress(compressed_array_cpu, compressed_size_step,
                decompressed_array_cpu, config, false);

            writer.Put<double>(variable_wt, (double *)decompressed_array_cpu, adios2::Mode::Sync);
            if (iterative_sim) {
                std::cout << "iterative update solution over time step\n";
                waveSim.init_u_np1((double *)decompressed_array_cpu);
            }
        } else {
            writer.Put<double>(variable_wt, waveSim.u_np1.data(), adios2::Mode::Sync);
        }
        iter_frame ++;
        writer.PerformPuts();
        writer.EndStep(); 
    }
    writer.Close();

    MPI_Finalize();
    return 0;
}













