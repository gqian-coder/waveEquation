#ifndef OBSTACLE_HPP
#define OBSTACLE_HPP

// multiple sphere objects
// radius: radius of the disk -- can be different for each
template <typename Real>
double* emulate_N_sphere(size_t Nx, size_t Ny, size_t Nz, size_t nDisks, size_t max_radius,
                        std::vector<size_t> x_pos, std::vector<size_t> y_pos, std::vector<size_t> z_pos);

template <typename Real>
double* emulate_two_slits(size_t Nx, size_t Ny, Real p1, Real p2,
                        Real width, Real dist_to_edge, Real thickness);

// two slits
template <typename Real>
double* emulate_two_slits(size_t Nx, size_t Ny);

#include "obstacle.tpp"
#endif
