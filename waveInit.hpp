#ifndef WAVEINIT_HPP
#define WAVEINIT_HPP

// A: magnitude of the sinusoid signal
// sigma: Gaussian sigma
// freq: frequency
template <typename Real>
void fun_gaussian_wave(Real *u, size_t Nx, size_t Ny, std::vector<Real> A,
                        std::vector<Real> sigma, std::vector<Real> freq_x,
                        std::vector<Real> freq_y, int nWaves);


// center_x, center_y: percentage, location in the x and y-axis
// A: magnitude of the sinusoid signal
// iter: iteration of the update
// freq: frequency -- 2*M_PI/T
template <typename Real>
void fun_sinusoidal(Real *u, size_t Nx, size_t Ny, float center_x, float center_y,
                    Real A, size_t iter, Real freq);


// A: magnitude
template <typename Real>
void fun_velocity(Real *u, size_t Nx, size_t Ny, Real dx, Real dy, Real A, Real sigma);


// initialization of a solid square
// NDx, NDy: square's dimension
template<typename Real>
void fun_square(Real *u, size_t Nx, size_t Ny, size_t NDx, size_t NDy, Real intensity);


// drop_probability: Raindrop probability (with each time tick) and intensity.
// NDx, NDy: droplet's region
// one droplet at each time
template <typename Real>
void fun_rainDrop(Real *u, size_t Nx, size_t Ny, size_t NDx, size_t NDy,
                  Real *gauss_template, float drop_probability);


// drop_probability: Raindrop probability (with each time tick) and intensity.
// drop multiple droplets
template <typename Real>
void fun_MultiRainDrop(Real *u, size_t Nx, size_t Ny, size_t NDx, size_t NDy,
                  Real *gauss_template, float drop_probability, int nDrops);

#include "waveInit.tpp"
#endif
