#include <cstdint>

#include <array>
#include <ostream>
#include <vector>

// gamma == 0 --> no dissipation
template<typename Real>
WaveEquation<Real>::WaveEquation(const std::size_t Nx, const std::size_t Ny, const std::size_t Nz,
                   Real dt, Real dh, Real c, Real gamma, int use_condition)
               : Nx(Nx), Ny(Ny), Nz(Nz), dt(dt), dh(dh), C(c), gamma(gamma), use_condition(use_condition)
{
    std::size_t data_size = (Nx+1) * (Ny+1) * (Nz+1);
    alpha  = C * dt / dh;
    alpha2 = alpha * alpha;
    u_n.resize(data_size, 0);
    u_nm1.resize(data_size, 0);
    u_np1.resize(data_size, 0);
}

template<typename Real>
void WaveEquation<Real>::init_u_n(Real *data)
{
    std::size_t data_size = (Nx+1) * (Ny+1) * (Nz+1);
    std::copy(data, data + data_size, u_n.begin());
}

template<typename Real>
void WaveEquation<Real>::init_u_np1(Real *data)
{
    std::size_t data_size = (Nx+1) * (Ny+1) * (Nz+1);
    std::copy(data, data + data_size, u_np1.begin());
}

template<typename Real>
void WaveEquation<Real>::update_3d()
{
    // old t -> new t-1
    std::copy(u_n.begin(), u_n.end(), u_nm1.begin());
    // old t+1 -> new t
    std::copy(u_np1.begin(), u_np1.end(), u_n.begin());
    // old t+1 <- 0
    std::fill(u_np1.begin(), u_np1.end(), 0);

    size_t offset_r, offset_c, k;
    Real gamma_t = gamma * dt;
    size_t dim1 = (Nz+1) * (Ny+1); 
    for (size_t r=1; r<Nx; r++) {
        offset_r = r * dim1;
        for (size_t c=1; c<Ny; c++) {
            offset_c = offset_r + c * (Nz+1); 
            for (size_t z=1; z<Nz; z++) {
                k = offset_c + z;
                u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2 * (
                           u_n[k-dim1] + u_n[k+dim1]
                         + u_n[k-Nz-1] + u_n[k+Nz+1] 
                         + u_n[k-1] + u_n[k+1] - 6*u_n[k]);
            }
        }
    }
    // Boundary conditions
    switch (use_condition) {
        case 1: // Dirichlet bound cond
        {
            // boundary @ r
            offset_r = Nx * dim1;
            for (size_t c=0; c<Ny+1; c++) {
                offset_c = c * (Nz+1);
                for (size_t z=0; z<Nz+1; z++) {
                    k = offset_c + z;
                    u_np1[k] = 0;
                    u_np1[k+offset_r] = 0;
                }
            }
            // boundary @ c
            offset_c = Ny * (Nz+1);
            for (size_t r=0; r<Nx+1; r++) {
                offset_r = r * dim1;
                for (size_t z=0; z<Nz+1; z++) {
                    u_np1[offset_r + z]    = 0;
                    u_np1[offset_r+ offset_c + z] = 0;
                }
            }
            // boundary @ z
            for (size_t r=0; r<Nx+1; r++) {
                offset_r = r * dim1;
                for (size_t c=0; c<Ny+1; c++) {
                    offset_c = offset_r + c * (Nz+1);
                    u_np1[offset_c]    = 0;
                    u_np1[offset_c+Nz] = 0;
                }
            }
            break;
        }
        case 2: // Mur boundary cond
        {
            Real kappa = (alpha - 1) / (alpha + 1);
            // boundary at R
            offset_r = Nx * dim1;
            for (size_t c=0; c<Ny+1; c++) {
                offset_c = c * (Nz+1);
                for (size_t z=0; z<Nz+1; z++) {
                    k  = offset_c + z;
                    u_np1[k] = u_n[k+dim1] + kappa * (u_np1[k+dim1] - u_n[k]);
                    k += offset_r; 
                    u_np1[k] = u_n[k-dim1] + kappa * (u_np1[k-dim1] - u_n[k]);
                }
            }
            // boundary at C
            for (size_t r=0; r<Nx+1; r++) {
                offset_r = r * dim1;
                for (size_t z=0; z<Nz+1; z++) {
                    k  = offset_r + z;  // c = 0
                    u_np1[k] = u_n[k+Nz+1] + kappa * (u_np1[k+Nz+1] - u_n[k]);
                    k += Ny * (Nz+1); // c = Ny
                    u_np1[k] = u_n[k-Nz-1] + kappa * (u_np1[k-Nz-1] - u_n[k]);
                }   
            }
            // boundary at Z
            for (size_t r=0; r<Nx+1; r++) {
                offset_r = r * dim1;
                for (size_t c=0; c<Ny+1; c++) {
                    k  = offset_r + c * (Nz+1); // z = 0
                    u_np1[k] = u_n[k+1] + kappa * (u_np1[k+1] - u_n[k]);
                    k += Nz;  // z = Nz
                    u_np1[k] = u_n[k-1] + kappa * (u_np1[k-1] - u_n[k]);
                }
            }
            break;
        }
        default:
            std::cout << "Invalid boundary condition option!" << std::endl;
    }
}

template<typename Real>
void WaveEquation<Real>::update_2d(bool init_cond)
{
    // old t -> new t-1
    std::copy(u_n.begin(), u_n.end(), u_nm1.begin());
    // old t+1 -> new t
    std::copy(u_np1.begin(), u_np1.end(), u_n.begin());
    // old t+1 <- 0
    std::fill(u_np1.begin(), u_np1.end(), 0);
    
    size_t offset_r, k;
    Real gamma_t = gamma * dt;
    for (size_t r=1; r<Nx; r++) {
        offset_r = r * (Ny+1);
        for (size_t c=1; c<Ny; c++) {
            k = offset_r + c;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2 * (
                       u_n[k-Ny-1] + u_n[k+Ny+1] + u_n[k-1] + u_n[k+1] - 4*u_n[k]);  
        }
    }
    
    // Boundary conditions
    switch (use_condition) {
        case 1: // Dirichlet bound cond
        {
            offset_r = Nx * (Ny+1);
            for (size_t c=0; c<Ny+1; c++) {
                u_np1[c] = 0;
                u_np1[c+offset_r] = 0;
            }
            for (size_t r=1; r<Nx; r++) {
                offset_r = r * (Ny+1);
                u_np1[offset_r]    = 0;
                u_np1[offset_r+Ny] = 0;
            }
            break;
        }
        case 2: // Mur boundary cond
        {
            Real kappa = (alpha - 1) / (alpha + 1); 
            offset_r = Nx * (Ny+1);
            for (size_t c=0; c<Ny+1; c++) {
                k = offset_r + c;
                u_np1[c] = u_n[Ny+1+c] + kappa * (u_np1[Ny+1+c] - u_n[c]);
                u_np1[k] = u_n[k-Ny-1] + kappa * (u_np1[k-Ny-1] - u_n[k]);
            }
            for (size_t r=0; r<Nx+1; r++) {
                k = r * (Ny+1);
                u_np1[k] = u_n[k+1] + kappa * (u_np1[k+1] - u_n[k]);
                k = k + Ny;
                u_np1[k] = u_n[k-1] + kappa * (u_np1[k-1] - u_n[k]);
            }
            break;
        }
        case 3: // Nuemann bound cond
        {
            // r,c = 0,0
            u_np1[0]  = u_n[0] + init_cond*(u_n[0]-u_nm1[0])  + 2*alpha2*(u_n[Ny+1] - 2*u_n[0] + u_n[1]);
            // r,c = 0,N_y
            u_np1[Ny] = u_n[Ny] + init_cond*(u_n[Ny]-u_nm1[Ny]) + 2*alpha2*(u_n[2*Ny+1] - 2*u_n[Ny] + u_n[Ny-1]);
            // r,c = Nx,0
            k = Nx * (Ny+1);
            u_np1[k]  = u_n[k] + init_cond*(u_n[k]-u_nm1[k]) + 2*alpha2*(u_n[k-Ny-1] - 2*u_n[k] + u_n[k+1]);
            // r,c = Nx,Ny
            k = offset_r + Ny;
            u_np1[k]  = u_n[k] + init_cond*(u_n[k]-u_nm1[k]) + 2*alpha2*(u_n[k-1] - 2*u_n[k] + u_n[k-1]);
            // r = 0
            for (size_t c=1; c<Ny-1; c++) {
                u_np1[c] = u_n[c] + init_cond*(u_n[c]-u_nm1[c]) + alpha2*(2*u_n[Ny+1+c] - 4*u_n[c] + u_n[c+1] + u_n[c-1]);
            }
            // c = 0
            for (size_t r=1; r<Nx-1; r++) {
                k = r * (Ny+1);
                u_np1[k] = u_n[k] + init_cond*(u_n[k]-u_nm1[k]) + alpha2*(u_n[k+Ny+1] - 4*u_n[k] + u_n[k-Ny-1] + 2*u_n[k+1]);
            }
            // r = Nx
            offset_r = Nx * (Ny+1);
            for (size_t c=1; c<Ny-1; c++) {
                k = offset_r + c;
                u_np1[k] = u_n[k] + init_cond*(u_n[k]-u_nm1[k]) + alpha2*(2*u_n[k-Ny-1] - 4*u_n[k] + u_n[k+1] + u_n[k-1]);
            }
            // c = Ny
            for (size_t r=1; r<Nx-1; r++) {
                offset_r = r * (Ny+1);
                k        = offset_r + Ny;
                u_np1[k] = u_n[k] + init_cond*(u_n[k]-u_nm1[k]) + alpha2*(u_n[k+Ny+1] - 4*u_n[k] + u_n[k-Ny-1] + 2*u_n[k-1]);
            }
            break;
        }
        case 4: // Periodic cond
        {
            // r = 0
            offset_r = (Nx-1) * (Ny+1);
            for (size_t c=1; c<Ny; c++) {
                u_np1[c] = 2 * u_n[c] - u_nm1[c] - gamma_t * (u_n[c]-u_nm1[c]) + alpha2 * (
                       u_n[offset_r+c] + u_n[c+Ny+1] + u_n[c-1] + u_n[c+1] - 4*u_n[c]);
            }
            // c = 0
            for (size_t r=1; r<Nx; r++) {
                k = r * (Ny+1);
                u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2 * (
                       u_n[k-Ny-1] + u_n[k+Ny+1] + u_n[k+Ny-1] + u_n[k+1] - 4*u_n[k]); 
            }
            // corner
            u_np1[0] = 2 * u_n[0] - u_nm1[0] - gamma_t * (u_n[0]-u_nm1[0]) + alpha2 * (
                       u_n[offset_r] + u_n[Ny+1] + u_n[Ny-1] + u_n[1] - 4*u_n[0]); 
            k = Nx*(Ny+1);
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2 * (
                       u_n[k-Ny-1] + u_n[Ny+1] + u_n[k+Ny-1] + u_n[k+1] - 4*u_n[k]);
            // open boundary
            offset_r = Nx * (Ny+1);
            for (size_t c=1; c<Ny; c++) {
                u_np1[c+offset_r] = u_np1[c] ;
            }
            for (size_t r=0; r<Nx+1; r++) {
                offset_r = r * (Ny+1);
                u_np1[offset_r+Ny] = u_np1[offset_r];
            }
            break;
        }
        default:
            std::cout << "Invalid boundary condition option!" << std::endl;
    }
}
