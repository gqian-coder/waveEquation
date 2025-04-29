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
#include <zfp.h>
#include "waveEquation.hpp"
#include "waveInit.hpp"
#include "obstacle.hpp"

unsigned char * zfp_compress(double * array, double tolerance, size_t r1, size_t r2, size_t r3, size_t *out_size){
    zfp_type type;     /* array scalar type */
    zfp_field* field;  /* array meta data */
    zfp_stream* zfp;   /* compressed stream */
    void* buffer;      /* storage for compressed stream */
    size_t bufsize;    /* byte size of compressed buffer */
    bitstream* stream; /* bit stream to write to or read from */
    size_t zfpsize;    /* byte size of compressed stream */

    /* allocate meta data for the 3D array a[nz][ny][nx] */
    type = zfp_type_double;
    if ((r1==1) && (r2==1)) {
        field = zfp_field_1d(array, type, r3);
    } else if (r1==1) {
        field = zfp_field_2d(array, type, r2, r3);
    }
    else {
        field = zfp_field_3d(array, type, r1, r2, r3);
    }
    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open(NULL);

    /* set compression mode and parameters via one of three functions */
    /*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
    /*  zfp_stream_set_precision(zfp, precision); */
    zfp_stream_set_accuracy(zfp, tolerance);

    /* allocate buffer for compressed data */
    bufsize = zfp_stream_maximum_size(zfp, field);
    buffer = malloc(bufsize);

    /* associate bit stream with allocated buffer */
    stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);

    zfpsize = zfp_compress(zfp, field);

    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
    *out_size = zfpsize;
    return (unsigned char *)buffer;
}

double * zfp_decompress(unsigned char * comp_data, double tolerance, size_t buffer_size, size_t r1, size_t r2, size_t r3){
    zfp_type type;     /* array scalar type */
    zfp_field* field;  /* array meta data */
    zfp_stream* zfp;   /* compressed stream */
    void* buffer;      /* storage for compressed stream */
    size_t bufsize;    /* byte size of compressed buffer */
    bitstream* stream; /* bit stream to write to or read from */

    /* allocate meta data for the 3D array a[nz][ny][nx] */
    double * array = (double *) malloc(r1 * r2 * r3 * sizeof(double));
    type = zfp_type_double;
    if ((r1==1) &&(r2==1)) {
        field = zfp_field_1d(array, type, r3);
    } else if (r1==1) {
        field = zfp_field_2d(array, type, r3, r2);
    }
    else
        field = zfp_field_3d(array, type, r3, r2, r1);

    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open(NULL);

    /* set compression mode and parameters via one of three functions */
    /*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
    /*  zfp_stream_set_precision(zfp, precision); */
    zfp_stream_set_accuracy(zfp, tolerance);

    /* allocate buffer for compressed data */
    bufsize = zfp_stream_maximum_size(zfp, field);
    // buffer = malloc(bufsize);
    buffer = (void *) comp_data;
    bufsize = buffer_size;

    /* associate bit stream with allocated buffer */
    stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);

    zfp_decompress(zfp, field);

    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
    return array;
}

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
    double Dz      = std::stof(argv[cnt_argv++]);
    // simulation spatial resolution
    double dh      = std::stof(argv[cnt_argv++]);

    size_t Nx = (size_t)std::ceil((double)Dx / dh);
    size_t Ny = (size_t)std::ceil((double)Dy / dh);
    size_t Nz = (size_t)std::ceil((double)Dz / dh);

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
        wave_c.resize(Nx*Ny*Nz);
        std::string v_fname(argv[cnt_argv++]);
        FILE *fp = fopen(v_fname.c_str(), "r");
        if (fp == NULL) {
            std::cout << "Error openning velocity mask file\n";
            return -1;
        }
        fread(wave_c.data(), sizeof(double), Nx*Ny*Nz, fp);
        fclose(fp);
    }
    double vp_max = *std::max_element(wave_c.begin(), wave_c.end());
    std::cout << "recommended timestep resolution < " << dh/vp_max/std::sqrt(2) << "\n";
    if (dt >= dh/vp_max/std::sqrt(2)) {
        std::cout << "Need smaller timestep resolution...\n";
        exit(-1);
    }
    // distribution of speed of sound: uniform | gaussian | mask
    std::string material_type(argv[cnt_argv++]);

    // dissipation rate
    double gamma   = std::stof(argv[cnt_argv++]); 
    float temp     = T / dt;
    size_t nframes = (size_t)(temp);
    // relative tolerance. tol = 0 means to write out uncompressed data after each timestep
    double tol     = std::stof(argv[cnt_argv++]);
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

    std::vector<size_t> dShape = {Nx, Ny, Nz};
    size_t cnt_data = Nx*Ny*Nz;

    std::cout << "simulating a domain of [" << Dx << "/" << Nx << ", " << Dy << "/" << Ny << ", " << Dz << "/" << Nz << "], at a spacing of " << dh << "\n";
    std::cout << "simulating " << nframes << " steps at a resolution of " << dt << "\n";
    std::cout << "initial function: ";
    if (init_fun==0) std::cout << "read from previous simulation steps in the file: " << fname << " at step " << init_ts << "\n";
    else if (init_fun==1) std::cout << "random raindrop\n";
    else if (init_fun==2) std::cout << "multiple rain drop\n";
    else if (init_fun==3) std::cout << "solid square\n";
    else if (init_fun==4) std::cout << "gaussian sinusoid waves\n";
    else if (init_fun==5) std::cout << "plane waves \n";
    
    if (obstacle_t==1) std::cout << "emulating multiple disk obstacles\n";
    if (obstacle_t==2) std::cout << "emulating two slits, recommending using plane wave\n";
    else std::cout << "no obstacle\n";
 
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
    } else if (!strcmp(material_type.c_str(), "mask")) {
        std::cout << "simulating material space using mask data\n"; 
    }
    std::cout << "tolerance = " << tol << " for u_n \n";

    if (src_on) std::cout << "turn on a  Gaussian derivative source with a relative ts = " << src_ts << "\n"; 
    
    adios2::ADIOS ad(MPI_COMM_WORLD);
    adios2::IO writer_io = ad.DeclareIO("Output");

    adios2::Engine writer = writer_io.Open(fname_wt, adios2::Mode::Write);
    adios2::Variable<double> variable_u;
    variable_u    = writer_io.DefineVariable<double>("u_data" , adios2::Dims{Nx, Ny, Nz}, adios2::Dims{0,0,0}, adios2::Dims{Nx, Ny, Nz});
    // store the non-zero masks
    size_t dim2 = Ny * Nz; 
    double max_intensity = 10.0;
    // for Gaussian wave
    double freq = 2.0 * M_PI / T;
    // for rain drops
    float drop_probability = 1;
    std::vector<double> gauss_template;
    double NDx = 12, NDy=24, NDz=24;
    size_t n_drops = 3;
    // for Gaussian pulse source
    double f0 = 100; // dominant frequency of the source (Hz)
    double t0 = 0.1; // source time shift (s) 
    size_t srcx = size_t(Dx * 0.5 / dh); 
    size_t srcy = size_t(Dy * 0.5 / dh);
    size_t srcz = size_t(Dz * 0.5 / dh);
    size_t src_pos = srcx * dim2 + srcy*Nz + srcz;
    double src_intensity = 50.0;

    size_t compressed_size_u = 0;
    double data_bytes = (double) Nx * Ny * Nz * sizeof(double);
    std::vector<double> compression_ratio(2);
    WaveEquation <double> waveSim(Nx-1, Ny-1, Nz-1, dt, dh, gamma, cfd_cond); 
    
    // initialize wave velocity
    std::vector<double> speed_sound(Nx*Ny*Nz);
    if (!strcmp(material_type.c_str(), "uniform")) {
        velocity_Layered_uniform(speed_sound.data(), wave_c, n_vp, Nx, Ny, Nz);
    } else if (!strcmp(material_type.c_str(), "gaussian")) {
        std::fill(speed_sound.begin(), speed_sound.end(), wave_c[0]);
        double peak_h = 0.5;
        velocity_Gaussian_3d(speed_sound.data(), wave_c[1], peak_h, Nx, Ny, Nz);
    } else if (!strcmp(material_type.c_str(), "sandwich")) {
        std::fill(speed_sound.begin(), speed_sound.end(), wave_c[0]);
        size_t width_slit = (size_t)(0.04 * (double)std::min(Ny, Nz));
        bool flg_v = velocity_sandwich(speed_sound.data(), wave_c[1], width_slit, Nx, Ny, Nz);
        if (flg_v==false) {
            std::cout << "Fail on creating a sandwich material space\n";
            exit(-1);
        }
    } else if (!strcmp(material_type.c_str(), "mask")) {
        velocity_mask(speed_sound.data(), wave_c.data(), Nx, Ny, Nz);
    }
    //FILE *fp = fopen("velocity_sandwich_3d.bin", "w");
    //fwrite(speed_sound.data(), sizeof(double), Nx*Ny*Nz, fp);
    //fclose(fp);
    std::cout << "velocity = " << *std::max_element(speed_sound.begin(), speed_sound.end()) << "\n";
    waveSim.init_vp(speed_sound.data());
    double *obstacle_m = NULL;
    if (obstacle_t) {
        switch (obstacle_t) {
            case 1: {
                size_t n_rods     = 10;
                size_t max_radius = size_t(0.015 * std::min(Ny, Nz));  
                size_t max_height = size_t(0.045 * Nx);
                obstacle_m        = emulate_N_rods<double>(Nx, Ny, Nz, n_rods, max_radius, max_height);
                break;
            }
            case 2: {
                size_t n_spheres  = 10;
                size_t max_radius = size_t(0.015 * std::min(Ny, Nz));
                obstacle_m        = emulate_N_sphere<double>(Nx, Ny, Nz, n_spheres, max_radius);
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
        NDz = (size_t) std::ceil(drop_width / dh);
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
        std::cout << "rainDrop region: " << NDx << " x " << NDy << " x" << NDz << " pixels\n";
        size_t radius_pow3 = NDx * NDy * NDz;
        for (size_t r=0; r<NDx; r++) {
            for (size_t c=0; c<NDy; c++) {
                for (size_t z=0; z<NDz; z++) {
                    if ((r-cx)*(r-cx)+(c-cy)*(c-cy)+(z-cz)*(z-cz) < radius_pow3) {
                        gauss_template[r*NDy*NDz+c*NDz+z] = max_intensity * exp(-(px[r]*px[r] + py[c]*py[c] + pz[z]*pz[z]));
                    }
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
            // load the checkpoint data compressed using non-optimized approach
            std::cout << "read u_{n}\n";
            variable = reader_io.InquireVariable<double>("u_data");
            
            std::cout << "total number of steps: " << variable.Steps() << ", read from " << init_ts << " timestep \n";
            std::vector<double> init_vpre(cnt_data);
            std::vector<double> init_v(cnt_data); 
            variable.SetSelection({adios2::Dims{0,0,0}, adios2::Dims{Nx, Ny, Nz}});
            variable.SetStepSelection({init_ts, 1}); 
            reader.Get(variable, init_vpre.data());
            reader.PerformGets();
            std::cout << "u_prev: min/max = " << *std::min_element(init_vpre.begin(), init_vpre.end()) << ", " << *std::max_element(init_vpre.begin(), init_vpre.end()) << "\n";
            std::cout << "begin to load the current step...\n";
            std::cout << "read u_{n+1}\n";
            variable.SetStepSelection({init_ts+1, 1});
            
            reader.Get(variable, init_v.data());
            reader.PerformGets();
            std::cout <<  "u_curr: min/max = " << *std::min_element(init_v.begin(), init_v.end()) << ", " << *std::max_element(init_v.begin(), init_v.end()) << "\n";
            
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
            }
            break;
        }
        case 1: { 
            fun_rainDrop<double>(waveSim.u_np1.data(), Nx, Ny, Nz, NDx, NDy, NDz, gauss_template.data(), drop_probability);
            break;
        }
        case 2: {
            fun_MultiRainDrop<double>(waveSim.u_np1.data(), Nx, Ny, Nz, NDx, NDy, NDz, gauss_template.data(), drop_probability, n_drops);
            break;
        }
        case 3: {
            fun_square<double>(waveSim.u_np1.data(), Nx, Ny, NDz, NDx, NDy, NDz, max_intensity);
            break;
        }
        case 4: {
            int n_waves = 4;
            std::vector<double> sigma = {50, 25, 35, 60};
            std::vector<double> intensity = {25, 30, 45, 15};
            std::vector<double> freq_x = {freq*T/18.0, freq*T/18.0, freq*T/24.0, freq*T/22.0};
            std::vector<double> freq_y = {freq*T/18.0, freq*T/32.0, freq*T/17.0, freq*T/18.0};
            std::vector<double> freq_z = {freq*T/18.0, freq*T/32.0, freq*T/17.0, freq*T/18.0};
            fun_gaussian_wave(waveSim.u_np1.data(), Nx, Ny, Nz, intensity, sigma, freq_x, freq_y, freq_z, n_waves);
            break;
        }
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
        std::cout << iter_frame << "/" << nframes << "\n";
        // source update
        if ((src_on) && (init_fun || iter_frame || ((tol==0) && (init_ts)>0))) {
            double ts = (init_fun==0) ? ((tol==0) ? 2 : 1) : 0; /* c.r. from t=2*/;
            ts += (double)(iter_frame + 1) + src_ts /* secondary c.r. */ + init_ts;
            ts  = ts * dt; 
            if (ts<0.25) { 
                waveSim.u_np1[src_pos] = src_Gaussian_pulse(f0, ts-t0, src_intensity);
                printf("update pos [%ld, %ld, %ld], src[%ld]=%.5e\n",srcx, srcy, srcz, (size_t)ts, waveSim.u_np1[src_pos]); 
            }
        }
        // wave update
        waveSim.update_3d();

        if (obstacle_t) {
            // apply Dirichlet boundary condition to the scattering on obstacles
            std::transform(waveSim.u_np1.begin(), waveSim.u_np1.end(), obstacle_m, waveSim.u_np1.begin(), std::multiplies<double>());
        }
        if ((iter_frame % wt_interval == 0) && (iter_frame>=cr_ts)) {
            writer.BeginStep();
            if (tol>0) {  // non-optimized compression
                unsigned char * compressed = zfp_compress(waveSim.u_n.data(), tol, Nx, Ny, Nz, &compressed_size_u); 
                std::cout << "compression ratios = " << data_bytes / (double) compressed_size_u << "\n";
                double *decompressed = zfp_decompress(compressed, tol, compressed_size_u, Nx, Ny, Nz); 
                writer.Put<double>(variable_u, (double *)decompressed, adios2::Mode::Sync);
                //double PE = potential_energy(waveSim.u_n.data(), decompressed, dShape);
                //printf("PE = %.8f, eb / sqrt(PE) = %.8f\n", std::sqrt(PE), tol / std::sqrt(PE));
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
    
    delete [] obstacle_m;
    MPI_Finalize();
    return 0;
}













