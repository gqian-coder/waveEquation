#include <random>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

// import velocity from a mask
template <typename Real>
void velocity_mask(Real *speed_sound, Real *mask, size_t Nx, size_t Ny, size_t Nz)
{
    size_t cnt_data = Nx*Ny*Nz;
    std::copy(mask, mask + cnt_data, speed_sound);
}

// wave velocity across domain, separating along x-axis, 1/n_vp fraction
template <typename Real>
void velocity_Layered_uniform(Real *speed_sound, std::vector<Real> wave_c, 
                                int n_vp, size_t Nx, size_t Ny, size_t Nz)
{
    std::vector<size_t> crossLayer(n_vp+1);
    for (int i=0; i<n_vp; i++) {
        crossLayer[i] = size_t((double)Nx / (double)n_vp * (double)i);
        // std::cout << i << ": " << crossLayer[i] << "\n";
    }
    crossLayer[n_vp] = Nx;
    size_t dim2 = Ny * Nz;
    size_t hdim, offset_r;
    for (int i=0; i<n_vp; i++) {
        hdim     = crossLayer[i+1] * dim2; 
        offset_r = crossLayer[i]   * dim2;
        // std::cout << hdim << ", " << offset_r << "\n";
        std::fill(speed_sound+offset_r, speed_sound+hdim, wave_c[i]); 
    }
} 

// Dual materials with the underneath layer in a bumped Gaussian shape 
template <typename Real> 
void velocity_Gaussian_2d(Real *speed_sound, Real wave_c, Real gaussian_peak, size_t Nx, size_t Ny)
{   
    Real mu = 0.0, sigma = 0.25;
    Real sigma2 = 2 * sigma * sigma;
    std::vector<Real> x_values(Ny), gaussian_values(Ny);
    Real step = 1.0 / (Real)(Ny-1);
    Real temp;
    gaussian_peak *= Nx;
    for (size_t i=0; i<Ny; i++) {
        x_values[i]        = -0.5 + i*step;
        temp = (x_values[i] - mu);
        gaussian_values[i] = gaussian_peak * std::exp(-(temp*temp) / sigma2); 
    }
    size_t j_end, offset_r;
    for (size_t i=0; i<Nx; i++) {
        for (size_t j=0; j<Ny; j++) {
            if (Nx-i<gaussian_values[j]) {
                //std::cout << "i = " << i << ", gv[" << j << "] = " << gaussian_values[j] << "\n";
                j_end = Ny-j+1; 
                offset_r = i * Ny;
                for (size_t k=j+offset_r; k<j_end+offset_r; k++) {
                    speed_sound[k] = wave_c;
                }
                break;
            }
        }
    }
}


// Dual materials with the underneath layer in a bumped Gaussian shape
template <typename Real>
void velocity_Gaussian_3d(Real *speed_sound, Real wave_c, Real gaussian_peak, 
                            size_t Nx, size_t Ny, size_t Nz)
{
    Real mu = 0.0, sigma = 0.25;
    Real sigma2 = 2 * sigma * sigma;
    std::vector<Real> x_values(Ny), y_values(Nz), gaussian_values(Nx);
    Real step_x = 1.0 / (Real)(Ny-1), step_y = 1.0 / (Real)(Nz-1);
    Real temp_x, temp_y;
    gaussian_peak *= Nx;
    size_t dim2 = Ny*Nz;
    size_t offset_r, offset_c;
    // generate gaussian values f(x,y) 
    for (size_t i=0; i<Ny; i++) {
        x_values[i] = -0.5 + i*step_x;
        offset_r = i * Nz; 
        temp_x = (x_values[i] - mu);
        temp_x = temp_x * temp_x; 
        for (size_t j=0; j<Nz; j++) {
            y_values[j] = -0.5 + j*step_y;
            temp_y = y_values[j] - mu;
            temp_y = temp_y * temp_y;
            gaussian_values[offset_r + j] = gaussian_peak * std::exp(-(temp_x + temp_y) / sigma2); 
        }
    }
    // filling the space under surface
    for (size_t i=0; i<Ny; i++) {
        offset_c = i*Nz;
        for (size_t j=0; j<Nz; j++) {
            for (size_t k=0; k<Nx; k++) {
                if (Nx-k<=gaussian_values[offset_c+j]) {
                    speed_sound[k*dim2+offset_c+j] = wave_c;
                }
            }
        }
    }
}


// sandwich materials with the middle layer different from the sides; divided by x-axis
template <typename Real>
bool velocity_sandwich(Real *speed_sound, Real wave_c, size_t width, 
                        size_t Nx, size_t Ny, size_t Nz)
{
    if (Ny <= width) {
        std::cout << "Need Ny>width to generate a sandwich shaped material space\n";    
        return false;
    } 
    size_t width_z = (Nz>width) ? width : 1;
    size_t dim1    = Ny * Nz;
    size_t offset_c, offset_r;
    size_t mid_c   = (size_t)(0.5 * (double)(Ny - width)); 
    size_t mid_z   = (Nz>width) ? (size_t)(0.5 * (double)(Nz - width)) : (Nz>2 ? (size_t)(0.5*(double)(Nz-1)):0);
    for (size_t r=0; r<Nx; r++) {
        offset_r = r *dim1;
        for (size_t c=mid_c; c<width+mid_c; c++) {
            offset_c = offset_r + c * Nz;
            for (size_t z=mid_z; z<width_z+mid_z; z++) {
                speed_sound[offset_c+z] = wave_c;
            }
        }
    }
    return true;
}

// Gaussian pulse point source: ts = time - t0
template <typename Real> 
Real src_Gaussian_pulse(Real freq, Real ts, Real A)  
{
    return (-2.0*A*ts * (freq*freq) * (std::exp(-(freq*freq) * (ts*ts)))) ;
}

// t0: nt * dt
template <typename Real>
void fun_Gaussian_pulse(Real *u, Real freq, Real t0, Real A, size_t xsrc, size_t ysrc, 
                    size_t zsrc, size_t Ny, size_t Nz)
{
    size_t k = xsrc * (Ny*Nz) + ysrc*Nz + zsrc;
    u[k] = src_Gaussian_pulse(freq, -t0, A); 
}


// an array of 2D plane wave starting from x=px 
template <typename Real>
void fun_plane_waves(Real *u, size_t Nx, size_t Ny, size_t px, 
                     Real A, Real freq, size_t nWaves)
{
    if (nWaves+px >= Nx) {
        std::cout << "Error::the number of waves cannot be larger than the size of data space\n";
    }
    double dy = 1.0 / (double)Ny;
    for (size_t r=px; r<nWaves; r++) {
        for (size_t c=0; c<Ny; c++) {
            u[r*Ny+c] = A * std::sin(1.0/nWaves*freq) * std::sin(M_PI*dy*c);
        }
    }
}


// Gaussian sinusoid waves
// A: magnitude of the sinusoid signal
// sigma: Gaussian sigma
// freq: frequency
template <typename Real>
void fun_gaussian_wave(Real *u, size_t Nx, size_t Ny, size_t Nz, std::vector<Real> A,
                        std::vector<Real> sigma, std::vector<Real> freq_x,
                        std::vector<Real> freq_y, std::vector<Real> freq_z, int nWaves)
{
    std::random_device rd;  // Non-deterministic random device (for seeding)
    std::mt19937 gen(rd()); // Mersenne Twister pseudo-random generator
    std::uniform_real_distribution<> dis(0.0, 1.0);

    size_t k, pos_x, pos_y, pos_z;
    float  center_x, center_y, center_z;
    std::vector <Real> px(Nx), py(Ny), pz(Nz);
    for (int i=0; i<nWaves; i++) {
        center_x = dis(gen);
        center_y = dis(gen);
        center_z = dis(gen);
        pos_x    = static_cast<size_t>(Nx*center_x);
        pos_y    = static_cast<size_t>(Ny*center_y);
        pos_z    = (Nz>1) ? static_cast<size_t>(Nz*center_z) : 0;
        std::cout << "pos_x = " << pos_x << ", pos_y = " << pos_y << ", pos_z = " << pos_z << "\n";
        std::cout << "freq_x = " << freq_x[i] << ", freq_y = " << freq_y[i] << ", freq_z = " << freq_z[i] << "\n";
        for (size_t r=0; r<Nx; r++) {
            px[r] = (static_cast<Real>(r)-pos_x);
            px[r] = px[r] * px[r];
        }
        for (size_t c=0; c<Ny; c++) {
            py[c] = (static_cast<Real>(c)-pos_y);
            py[c] = py[c] * py[c];
        }
        for (size_t z=0; z<Nz; z++) {
            pz[z] = (static_cast<Real>(z)-pos_z);
            pz[z] = pz[z] * pz[z];
        }
        Real sigma_2 = 2*sigma[i]*sigma[i];
        size_t dim1  = Ny * Nz;
        size_t offset_x, offset_y;
        for (size_t r=0; r<Nx; r++) {
            offset_x = r * dim1;
            for (size_t c=0; c<Ny; c++) {
                offset_y = c * Nz;
                for (size_t z=0; z<Nz; z++) {
                    k = offset_x + offset_y + z;
                    u[k] += A[i] * exp(-(px[r] + py[c] + pz[z]) / sigma_2) * sin(freq_x[i]*r + freq_y[i]*c +freq_z[i]*z);
                }
            }
        }
    }
}


// A: magnitude
template <typename Real>
void fun_Gaussian(Real *u, size_t Nx, size_t Ny, size_t Nz,
                  Real dx, Real dy, Real dz, Real A, Real sigma)
{
    size_t offset_x, offset_y, k;
    std::vector <Real> px(Nx), py(Ny), pz(Nz);
    for (size_t r=0; r<Nx; r++) {
        px[r] = (static_cast<Real>(r)-Nx)*dx;
        px[r] = px[r] * px[r];
    }
    for (size_t c=0; c<Ny; c++) {
        py[c] = (static_cast<Real>(c)-Ny)*dy;
        py[c] = py[c] * py[c];
    }
    for (size_t z=0; z<Nz; z++) {
        pz[z] = Nz>1 ? (static_cast<Real>(z)-Nz)*dz : 0;
        pz[z] = pz[z] * pz[z];
    }
    Real sigma_2 = 2*sigma*sigma;
    size_t dim1  = Ny * Nz;
    for (size_t r=0; r<Nx; r++) {
        offset_x = r * dim1;
        for (size_t c=0; c<Ny; c++) {
            offset_y = c * Nz;
            for (size_t z=0; z<Nz; z++) {
                k  = offset_x + offset_y + z;
                u[k] = A * exp(-((px[r] + py[c] + pz[z]) / sigma_2));
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
    float random_x = 0.5;//dis(gen);
    float random_y = 0.5;//dis(gen);
    float random_z = 0.5;//dis(gen);
    size_t x = static_cast<size_t>(random_x * (Nx-cx-1) + cx);
    size_t y = static_cast<size_t>(random_y * (Ny-cy-1) + cy);
    size_t z = static_cast<size_t>(random_z * (Nz-cz-1) + cz);
    size_t dim1  = Ny * Nz;
    size_t offset_x, offset_y, k;
    for (size_t r=0; r<NDx; r++) {
        offset_x = (x-cx+r)*dim1;
        for (size_t c=0; c<NDy; c++) {
            offset_y = (y-cy+c)*Nz;
            for (size_t h=0; h<NDz; h++) {
                k = offset_x + offset_y + z-cz+h;
                u[k] = intensity;
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
    std::cout << cx << ", " << cy << ", " << cz << "\n";
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
        // std::cout << "x, y, z = " << x << ", " << y  << ", "<< cx << ", " << cy << ", " << cz << ", " << NDx << ", " << NDy << ", " << NDz << Nx << ", " << Ny<< ", " << Nz << "\n";
        size_t dim1  = Ny * Nz;
        size_t dimD1 = NDy * NDz;
        size_t offset_x, offset_y, k;
        for (size_t r=0; r<NDx; r++) {
            offset_x = (x-cx+r)*dim1;
            for (size_t c=0; c<NDy; c++) {
                offset_y = (y-cy+c)*Nz;
                for (size_t h=0; h<NDz; h++) {
                    k = offset_x + offset_y + z-cz+h;
                    u[k] += gauss_template[r*dimD1+c*NDz+h];
                }
            }
        }
    }
}

// drop_probability: Raindrop probability (with each time tick) and intensity.
// drop multiple droplets
template <typename Real>
void fun_MultiRainDrop(Real *u, size_t Nx, size_t Ny, size_t Nz, 
                       size_t NDx, size_t NDy, size_t NDz,
                       Real *gauss_template, float drop_probability, size_t nDrops)
{
    size_t cx = (size_t)(NDx/2);
    size_t cy = (size_t)(NDy/2);
    size_t cz = (size_t)(NDz/2);
    std::random_device rd;  // Non-deterministic random device (for seeding)
    std::mt19937 gen(rd()); // Mersenne Twister pseudo-random generator
    double noRainZone = 0.1;
    std::uniform_real_distribution<> dis(0.0+noRainZone, 1.0-noRainZone);
    float random_number = dis(gen);
    //std::cout << "probability: " << random_number << "\n";
    // try not to generate rain drops closing to edges
    std::vector<float> px  = {0.275, 0.65, 0.35, 0.7, 0.8};
    std::vector<float> py  = {0.35, 0.3, 0.725, 0.65, 0.3};
    std::vector<float> pz  = {0.34, 0.56, 0.75, 0.25, 0.58};
    std::vector<float> mag = {0.7, 0.5, 0.85, 0.6, 0.9};
    for (size_t d=0; d<nDrops; d++) {
        if (random_number < drop_probability) {
            float random_x = px[d];//dis(gen);
            float random_y = py[d];//dis(gen);
            float random_z = pz[d];//dis(gen);
            size_t x = static_cast<size_t>(random_x * (Nx-cx-1) + cx);
            size_t y = static_cast<size_t>(random_y * (Ny-cy-1) + cy);
            size_t z = static_cast<size_t>(random_z * (Nz-cz-1) + cz);
            std::cout << "x, y, z = " << x << ", " << y  << ", "<< z << ", " << cx << ", " << cy << ", " << cz << ", " << NDx << ", " << NDy << ", " << NDz << ", " << Nx << ", " << Ny<< ", " << Nz << "\n";
            size_t dim1  = Ny * Nz;
            size_t dimD1 = NDy * NDz;
            size_t offset_x, offset_y, k;
            float intensity = mag[d];//dis(gen);
            for (size_t r=0; r<NDx; r++) {
                offset_x = (x-cx+r)*dim1;
                for (size_t c=0; c<NDy; c++) {
                    offset_y = (y-cy+c)*Nz;
                    for (size_t h=0; h<NDz; h++) {
                        k = offset_x + offset_y + z-cz+h;
                        u[k] += intensity*gauss_template[r*dimD1+c*NDz+h];
                        //std::cout << gauss_template[r*NDy+c] << "\n";
                    }
                }
            }
        }
    }
    //std::cout << *std::min_element(u, u+Nx*Ny*Nz) << ", " << *std::max_element(u, u+Nx*Ny*Nz) << "\n";
}

