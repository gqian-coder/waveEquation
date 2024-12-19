// drop_probability: Raindrop probability (with each time tick) and intensity.
#include <random>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

// multiple disk objects
// radius: radius of the disk -- can be different for each
template <typename Real>
double* emulate_N_disk(size_t Nx, size_t Ny, size_t nDisks, size_t max_radius,
                        std::vector<size_t> x_pos, std::vector<size_t> y_pos)
{
    std::random_device rd;  // Non-deterministic random device (for seeding)
    std::mt19937 gen(rd()); // Mersenne Twister pseudo-random generator
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    std::vector<int> radius(nDisks);    
    double *mask = new double [Nx*Ny];
    std::fill(mask, mask+Nx*Ny, 1.0);
    
    std::vector<size_t> x_prev(nDisks, Nx), y_prev(nDisks, Ny);
    size_t x, y;
    for (size_t d=0; d<nDisks; d++) {
        //radius[d]      = size_t(dis(gen) * max_radius);
        radius[d] = max_radius;
        if ((x_pos.size()<nDisks) || (y_pos.size()<nDisks)) { 
            float random_x = dis(gen);
            float random_y = dis(gen);
            x = static_cast<size_t>(random_x * (Nx-radius[d]-1) + radius[d]);
            y = static_cast<size_t>(random_y * (Ny-radius[d]-1) + radius[d]);
            // check for overlapping
            bool nonOverlap = true;
            for (size_t m=0; m<d; m++) { 
                nonOverlap &= ((std::abs((double)x-(double)x_prev[m]) > radius[d]+radius[m]) || (std::abs((double)y-(double)y_prev[m]) > radius[d]+radius[m])); 
            }
            while (nonOverlap==false) {
                nonOverlap = true;
                random_x = dis(gen);
                random_y = dis(gen);
                x = static_cast<size_t>(random_x * (Nx-radius[d]-1) + radius[d]);
                y = static_cast<size_t>(random_y * (Ny-radius[d]-1) + radius[d]);
                for (size_t m=0; m<d; m++) {
                    nonOverlap = nonOverlap & ((std::abs((double)x-(double)x_prev[m]) > radius[d]+radius[m]) || (std::abs((double)y-(double)y_prev[m]) < radius[d]+radius[m]));
                }
            }
            x_prev[d] = x; 
            y_prev[d] = y; 
        } else {
            x = x_pos[d];
            y = y_pos[d];
        }
        std::cout << "obstacle " << d << ": radius = " << radius[d] << ", pos {" << x << ", " << y << "}\n"; 
        int diam = radius[d] * 2;
        int radius_sqrt = radius[d] * radius[d];
        for (int r=0; r<diam; r++) {
            for (int c=0; c<diam; c++) {
                if ((r-radius[d])*(r-radius[d]) + (c-radius[d])*(c-radius[d]) < radius_sqrt) {
                    mask[(x-radius[d]+r)*Ny + (y-radius[d]+c)] = 0.0; 
                }
            }
        }
    }
    std::cout << "finish obstacle construction\n";
    return mask;
}

 
// two slits
template <typename Real>
double* emulate_two_slits(size_t Nx, size_t Ny)
{

    double *mask = new double[Nx*Ny];
    std::fill(mask, mask+Nx*Ny, 1.0);

    // wall
    size_t thickness    = size_t ((1.0/32.0) * Nx);
    size_t dist_to_edge = size_t ((1.0/5.0) * Nx);
    for (size_t r=dist_to_edge; r<dist_to_edge+thickness; r++) {
        for (size_t c=0; c<Ny; c++) {
            mask[r*Ny+c] = 0.0;
        }
    } 
    // slits
    size_t width = size_t ((1.0/16.0) * Ny);
    size_t p1     = size_t ((5.0/16.0) * Ny);
    size_t p2     = size_t ((5.0/8.0 ) * Ny);
    for (size_t r=dist_to_edge; r<dist_to_edge+thickness; r++) {
        for (size_t c=p1; c<p1+width; c++) {
            mask[r*Ny+c] = 1.0;
        }
        for (size_t c=p2; c<p2+width; c++) {
            mask[r*Ny+c] = 1.0;
        }
    }
    return mask;
}

