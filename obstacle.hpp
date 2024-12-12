#ifndef OBSTACLE_HPP
#define OBSTACLE_HPP

// multiple disk objects
// radius: radius of the disk -- can be different for each
template <typename Real>
double* emulate_N_disk(size_t Nx, size_t Ny, size_t nDisks, size_t radius);

// two slits
template <typename Real>
double* emulate_two_slits(size_t Nx, size_t Ny);

#include "obstacle.tpp"
#endif
