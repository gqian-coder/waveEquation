#include <random>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

// an array of plane wave starting from x=0 
template <typename Real>
void fun_plane_waves(Real *u, size_t Nx, size_t Ny, Real A, Real freq, size_t nWaves)
{
    if (nWaves >= Nx) {
        std::cout << "Error::the number of waves cannot be larger than the size of data space\n";
    }
    double dy = 1.0 / (double)Ny;
    for (size_t r=0; r<nWaves; r++) {
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
void fun_gaussian_wave(Real *u, size_t Nx, size_t Ny, std::vector<Real> A,
                        std::vector<Real> sigma, std::vector<Real> freq_x,
                        std::vector<Real> freq_y, int nWaves)
{
    std::random_device rd;  // Non-deterministic random device (for seeding)
    std::mt19937 gen(rd()); // Mersenne Twister pseudo-random generator
    std::uniform_real_distribution<> dis(0.0, 1.0);

    size_t offset_x, k, pos_x, pos_y;
    float  center_x, center_y;
    std::vector <Real> px(Nx), py(Ny);
    for (int i=0; i<nWaves; i++) {
        center_x = dis(gen);
        center_y = dis(gen);
        pos_x    = static_cast<size_t>(Nx*center_x);
        pos_y    = static_cast<size_t>(Ny*center_y);
        std::cout << "pos_x = " << pos_x << ", pos_y = " << pos_x << "\n";
        std::cout << "freq_x = " << freq_x[i] << ", freq_y = " << freq_y[i] << "\n";
        for (size_t r=0; r<Nx; r++) {
            px[r] = (static_cast<Real>(r)-pos_x);
            px[r] = px[r] * px[r];
        }
        for (size_t c=0; c<Ny; c++) {
            py[c] = (static_cast<Real>(c)-pos_y);
            py[c] = py[c] * py[c];
        }
        Real sigma_2 = 2*sigma[i]*sigma[i];
        for (size_t r=0; r<Nx; r++) {
            offset_x = r * Ny;
            for (size_t c=0; c<Ny; c++) {
                k = offset_x + c;
                u[k] += A[i] * exp(-(px[r] + py[c]) / sigma_2) * sin(freq_x[i]*r + freq_y[i]*c);
            }
        }
    }
}

// center_x, center_y: percentage, location in the x and y-axis
// A: magnitude of the sinusoid signal
// iter: iteration of the update
// freq: frequency -- 2*M_PI/T
template <typename Real>
void fun_sinusoidal(Real *u, size_t Nx, size_t Ny, float center_x, float center_y,
                    Real A, size_t iter, Real freq)
{
    if ((center_x <= 0) || (center_x>=1) || (center_y<=0) || (center_y>=1)) {
        std::cout << "The initial position of the sinusoidal ring must be a positive number and less than the boundary\n";
        return;
    }
    size_t pos_x = static_cast<size_t>(Nx*center_x);
    size_t pos_y = static_cast<size_t>(Ny*center_y);
    u[pos_x * Ny + pos_y] += A * sin(iter*freq);
    //std::cout << pos_x * Ny + pos_y << ", " << pos_x << ", " << Ny << ", " << pos_y << ", " << u[pos_x * Ny + pos_y] << "\n";
}

// A: magnitude
template <typename Real>
void fun_velocity(Real *u, size_t Nx, size_t Ny, Real dx, Real dy, Real A, Real sigma)
{
    size_t offset_x, k;
    std::vector <Real> px(Nx), py(Ny);
    for (size_t r=0; r<Nx; r++) {
        px[r] = (static_cast<Real>(r)-Nx)*dx;
    }
    for (size_t c=0; c<Ny; c++) {
        py[c] = (static_cast<Real>(c)-Ny)*dy;
    }
    for (size_t r=0; r<Nx; r++) {
        offset_x = r * Ny;
        for (size_t c=0; c<Ny; c++) {
            k  = offset_x + c;
            u[k] = A * exp(-((px[r]*px[r] + py[c]*py[c])/sigma));
        }
    }
}

// initialization of a solid square
// NDx, NDy: square's dimension
template<typename Real>
void fun_square(Real *u, size_t Nx, size_t Ny, size_t NDx, size_t NDy, Real intensity)
{
    std::cout << "fun_square\n";
    size_t cx = (size_t)(NDx/2);
    size_t cy = (size_t)(NDy/2);
    std::random_device rd;  // Non-deterministic random device (for seeding)
    std::mt19937 gen(rd()); // Mersenne Twister pseudo-random generator
    std::uniform_real_distribution<> dis(0.0, 1.0);
    float random_x = 0.5;//dis(gen);
    float random_y = 0.5;//dis(gen);
    size_t x = static_cast<size_t>(random_x * (Nx-cx-1) + cx);
    size_t y = static_cast<size_t>(random_y * (Ny-cy-1) + cy);
    for (size_t r=0; r<NDx; r++) {
        for (size_t c=0; c<NDy; c++) {
            u[(x-cx + r)*Ny + (y-cy+c)] = intensity;
            //std::cout << gauss_template[r*NDy+c] << "\n";
        }
    }
}

// drop_probability: Raindrop probability (with each time tick) and intensity.
// NDx, NDy: droplet's region
// one droplet at each time
template <typename Real>
void fun_rainDrop(Real *u, size_t Nx, size_t Ny, size_t NDx, size_t NDy,
                  Real *gauss_template, float drop_probability)
{
    size_t cx = (size_t)(NDx/2);
    size_t cy = (size_t)(NDy/2);
    std::random_device rd;  // Non-deterministic random device (for seeding)
    std::mt19937 gen(rd()); // Mersenne Twister pseudo-random generator
    std::uniform_real_distribution<> dis(0.0, 1.0);
    float random_number = dis(gen);
    //std::cout << "probability: " << random_number << "\n";
    if (random_number < drop_probability) {
        float random_x = 0.5;//dis(gen);
        float random_y = 0.5;//dis(gen);
        size_t x = static_cast<size_t>(random_x * (Nx-cx-1) + cx);
        size_t y = static_cast<size_t>(random_y * (Ny-cy-1) + cy);
        //std::cout << "x, y = " << x << ", " << y  << ", "<< cx << ", " << cy << ", " << NDx << ", " << NDy << ", " << Nx << ", " << Ny<< "\n";
        for (size_t r=0; r<NDx; r++) {
            for (size_t c=0; c<NDy; c++) {
                u[(x-cx + r)*Ny + (y-cy+c)] += gauss_template[r*NDy+c];
                //std::cout << gauss_template[r*NDy+c] << "\n";
            }
        }
    }
}

// drop_probability: Raindrop probability (with each time tick) and intensity.
// drop multiple droplets
template <typename Real>
void fun_MultiRainDrop(Real *u, size_t Nx, size_t Ny, size_t NDx, size_t NDy,
                  Real *gauss_template, float drop_probability, size_t nDrops)
{
    size_t cx = (size_t)(NDx/2);
    size_t cy = (size_t)(NDy/2);
    std::random_device rd;  // Non-deterministic random device (for seeding)
    std::mt19937 gen(rd()); // Mersenne Twister pseudo-random generator
    std::uniform_real_distribution<> dis(0.0, 1.0);
    float random_number = dis(gen);
    //std::cout << "probability: " << random_number << "\n";
    for (size_t d=0; d<nDrops; d++) {
        if (random_number < drop_probability) {
            float random_x = dis(gen);
            float random_y = dis(gen);
            size_t x = static_cast<size_t>(random_x * (Nx-cx-1) + cx);
            size_t y = static_cast<size_t>(random_y * (Ny-cy-1) + cy);
            std::cout << "x, y = " << x << ", " << y  << ", "<< cx << ", " << cy << ", " << NDx << ", " << NDy << ", " << Nx << ", " << Ny<< "\n";
            for (size_t r=0; r<NDx; r++) {
                for (size_t c=0; c<NDy; c++) {
                    float intensity = dis(gen);
                    u[(x-cx + r)*Ny + (y-cy+c)] += intensity*gauss_template[r*NDy+c];
                    //std::cout << gauss_template[r*NDy+c] << "\n";
                }
            }
        }
    }
    //std::cout << *std::min_element(u, u+Nx*Ny) << ", " << *std::max_element(u, u+Nx*Ny) << "\n";
}

