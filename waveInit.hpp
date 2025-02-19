#ifndef WAVEINIT_HPP
#define WAVEINIT_HPP

// an array of plane wave starting from x=px
template <typename Real>
void fun_plane_waves(Real *u, size_t Nx, size_t Ny, size_t px,
                     Real A, Real freq, size_t nWaves);

// A: magnitude of the sinusoid signal
// sigma: Gaussian sigma
// freq: frequency
template <typename Real>
void fun_gaussian_wave(Real *u, size_t Nx, size_t Ny, size_t Nz, std::vector<Real> A,
                        std::vector<Real> sigma, std::vector<Real> freq_x,
                        std::vector<Real> freq_y, std::vector<Real> freq_z, int nWaves);


// A: magnitude
template <typename Real>
void fun_Gaussian(Real *u, size_t Nx, size_t Ny, size_t Nz,
                  Real dx, Real dy, Real dz, Real A, Real sigma);


// initialization of a solid square
// NDx, NDy: square's dimension
template<typename Real>
void fun_square(Real *u, size_t Nx, size_t Ny, size_t Nz,
                size_t NDx, size_t NDy, size_t NDz, Real intensity);

// drop_probability: Raindrop probability (with each time tick) and intensity.
// NDx, NDy: droplet's region
// one droplet at each time
template <typename Real>
void fun_rainDrop(Real *u, size_t Nx, size_t Ny, size_t Nz,
                  size_t NDx, size_t NDy, size_t NDz,
                  Real *gauss_template, float drop_probability);

// drop_probability: Raindrop probability (with each time tick) and intensity.
// drop multiple droplets
template <typename Real>
void fun_MultiRainDrop(Real *u, size_t Nx, size_t Ny, size_t Nz,
                       size_t NDx, size_t NDy, size_t NDz,
                       Real *gauss_template, float drop_probability, size_t nDrops);

// Gaussian pulse point source: ts = time - t0
template <typename Real>
Real src_Gaussian_pulse(Real freq, Real ts, Real A);

// initializing the domain using a Gaussian pulse
template <typename Real>
void fun_Gaussian_pulse(Real *u, Real freq, Real t0, Real A, size_t xsrc, size_t ysrc,
                    size_t zsrc, size_t Ny, size_t Nz);

// wave velocity across domain, separating along x-axis, 1/n_vp fraction
template <typename Real>
void velocity_Layered_uniform(Real *speed_sound, std::vector<Real> wave_c,
                                int n_vp, size_t Nx, size_t Ny, size_t Nz);

template <typename Real>
void velocity_Gaussian_2d(Real *speed_sound, Real wave_c, Real gaussian_peak, size_t Nx, size_t Ny);

template <typename Real>
void velocity_Gaussian_3d(Real *speed_sound, Real wave_c, Real gaussian_peak,
                            size_t Nx, size_t Ny, size_t Nz);

template <typename Real>
bool velocity_sandwich(Real *speed_sound, Real wave_c, size_t width,
                        size_t Nx, size_t Ny, size_t Nz);

// import velocity from a mask
template <typename Real>
void velocity_mask(Real *speed_sound, Real *mask, size_t Nx, size_t Ny, size_t Nz);

#include "waveInit.tpp"
#endif
