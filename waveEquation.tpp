#include <cstdint>

#include <array>
#include <ostream>
#include <vector>

// gamma == 0 --> no dissipation
template<typename Real>
WaveEquation<Real>::WaveEquation(const std::size_t Nx, const std::size_t Ny, const std::size_t Nz,
                   Real dt, Real dh, Real gamma, int use_condition)
               : Nx(Nx), Ny(Ny), Nz(Nz), dt(dt), dh(dh), gamma(gamma), use_condition(use_condition)
{
    std::size_t data_size = (Nx+1) * (Ny+1) * (Nz+1);
    u_n.resize(data_size, 0);
    u_nm1.resize(data_size, 0);
    u_np1.resize(data_size, 0);
    vp.resize(data_size, 0);
    alpha2.resize(data_size, 0);
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
void WaveEquation<Real>::init_vp(Real *data)
{
    std::size_t data_size = (Nx+1) * (Ny+1) * (Nz+1);
    std::copy(data, data + data_size, vp.begin());

    Real delta = dt / dh ;
    for (size_t i=0; i<data_size; i++) {
        alpha2[i] = vp[i] * delta;
        alpha2[i] = alpha2[i] * alpha2[i]; 
    }
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
                u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
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
            Real delta = dt / dh;
            Real kappa;
            // boundary at R
            offset_r = Nx * dim1;
            for (size_t c=0; c<Ny+1; c++) {
                offset_c = c * (Nz+1);
                for (size_t z=0; z<Nz+1; z++) {
                    k  = offset_c + z;
                    kappa    = (vp[k]*delta - 1) / (vp[k]*delta + 1);
                    u_np1[k] = u_n[k+dim1] + kappa * (u_np1[k+dim1] - u_n[k]);
                    k += offset_r; 
                    kappa    = (vp[k]*delta - 1) / (vp[k]*delta + 1);
                    u_np1[k] = u_n[k-dim1] + kappa * (u_np1[k-dim1] - u_n[k]);
                }
            }
            // boundary at C
            for (size_t r=0; r<Nx+1; r++) {
                offset_r = r * dim1;
                for (size_t z=0; z<Nz+1; z++) {
                    k  = offset_r + z;  // c = 0
                    kappa    = (vp[k]*delta - 1) / (vp[k]*delta + 1);
                    u_np1[k] = u_n[k+Nz+1] + kappa * (u_np1[k+Nz+1] - u_n[k]);
                    k += Ny * (Nz+1); // c = Ny
                    kappa    = (vp[k]*delta - 1) / (vp[k]*delta + 1);
                    u_np1[k] = u_n[k-Nz-1] + kappa * (u_np1[k-Nz-1] - u_n[k]);
                }   
            }
            // boundary at Z
            for (size_t r=0; r<Nx+1; r++) {
                offset_r = r * dim1;
                for (size_t c=0; c<Ny+1; c++) {
                    k  = offset_r + c * (Nz+1); // z = 0
                    kappa    = (vp[k]*delta - 1) / (vp[k]*delta + 1);
                    u_np1[k] = u_n[k+1] + kappa * (u_np1[k+1] - u_n[k]);
                    k += Nz;  // z = Nz
                    kappa    = (vp[k]*delta - 1) / (vp[k]*delta + 1);
                    u_np1[k] = u_n[k-1] + kappa * (u_np1[k-1] - u_n[k]);
                }
            }
            break;
        }
        case 3: // Periodic cond
        // r = 0
        offset_r = Nx * dim1;
        for (size_t c=1; c<Ny; c++) {
            offset_c = c * (Nz+1);
            for (size_t z=1; z<Nz; z++) {
                k = offset_c + z; 
                u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[offset_r+k] + u_n[k+dim1]
                         + u_n[k-Nz-1] + u_n[k+Nz+1]
                         + u_n[k-1] + u_n[k+1] - 6*u_n[k]);
            }
        }
        // c = 0
        offset_c = Ny * (Nz+1);
        for (size_t r=1; r<Nx; r++) {
            offset_r = r * dim1;
            for (size_t z=1; z<Nz; z++) {
                k = offset_r + z;
                u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[k+dim1]
                         + u_n[k+offset_c] + u_n[k+Nz+1]
                         + u_n[k-1] + u_n[k+1] - 6*u_n[k]);   
            }
        }
        // z = 0
        for (size_t r=1; r<Nx; r++) {
            offset_r = r * dim1;
            for (size_t c=1; c<Ny; c++) {
                k = offset_r + c*(Nz+1);
                u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[k+dim1]
                         + u_n[k-Nz-1] + u_n[k+Nz+1]
                         + u_n[k+Nz] + u_n[k+1] - 6*u_n[k]);
            }
        }
        // r = Nx
        offset_r = Nx * dim1;
        for (size_t c=1; c<Ny; c++) {
            offset_c = c * (Nz+1);
            for (size_t z=1; z<Nz; z++) {
                k = offset_r + offset_c + z;
                u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[offset_c+z]
                         + u_n[k-Nz-1] + u_n[k+Nz+1]
                         + u_n[k-1] + u_n[k+1] - 6*u_n[k]);
            }
        }
        // c = Ny
        offset_c = Ny * (Nz+1);
        for (size_t r=1; r<Nx; r++) {
            offset_r = r * dim1;
            for (size_t z=1; z<Nz; z++) {
                k = offset_r + offset_c + z;
                u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[k+dim1]
                         + u_n[k-Nz-1] + u_n[offset_r+z]
                         + u_n[k-1] + u_n[k+1] - 6*u_n[k]);
            }
        }
        // z = Nz
        for (size_t r=1; r<Nx; r++) {
            offset_r = r * dim1;
            for (size_t c=1; c<Ny; c++) {
                k = offset_r + c*(Nz+1) + Nz;
                u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[k+dim1]
                         + u_n[k-Nz-1] + u_n[k+Nz+1]
                         + u_n[k-1] + u_n[k-Nz] - 6*u_n[k]);
            }
        }
        // r = 0, c = 0 or Ny
        offset_c = Ny * (Nz+1);
        offset_r = Nx * dim1;
        for (size_t z=1; z<Nz; z++) {
            k = z;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k+offset_r] + u_n[k+dim1]
                         + u_n[k+offset_c] + u_n[k+Nz+1]
                         + u_n[k-1] + u_n[k+1] - 6*u_n[k]);
            k += offset_c;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k+offset_r] + u_n[k+dim1]
                         + u_n[k-Nz-1] + u_n[z]
                         + u_n[k-1] + u_n[k+1] - 6*u_n[k]);   
        }
        // r = Nx, c = 0 or Ny
        for (size_t z=1; z<Nz; z++) {
            k = offset_r + z;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[z]
                         + u_n[k+offset_c] + u_n[k+Nz+1]
                         + u_n[k-1] + u_n[k+1] - 6*u_n[k]);
            k += offset_c;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[z]
                         + u_n[k-Nz-1] + u_n[k-offset_c]
                         + u_n[k-1] + u_n[k+1] - 6*u_n[k]);
        } 
        // c = 0, z = 0 or Nz
        for (size_t r=1; r<Nx; r++) {
            k = r * dim1;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[k+dim1]
                         + u_n[k+offset_c] + u_n[k+Nz+1]
                         + u_n[k+Nz] + u_n[k+1] - 6*u_n[k]);
            k += Nz; 
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[k+dim1]
                         + u_n[k+offset_c] + u_n[k+Nz+1]
                         + u_n[k-1] + u_n[k-Nz] - 6*u_n[k]);
        } 
        // c = Ny, z = 0 or Nz
        for (size_t r=1; r<Nx; r++) {
            k = r * dim1 + offset_c;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[k+dim1]
                         + u_n[k-Nz-1] + u_n[k-offset_c]
                         + u_n[k+Nz] + u_n[k+1] - 6*u_n[k]);
            k += Nz;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[k+dim1]
                         + u_n[k-Nz-1] + u_n[k-offset_c]
                         + u_n[k-1] + u_n[k-Nz] - 6*u_n[k]);
        }
        // z = 0, r = 0 or Nx
        for (size_t c=1; c<Ny; c++) {
            k = c * (Ny+1);
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k+offset_r] + u_n[k+dim1]
                         + u_n[k-Nz-1] + u_n[k+Nz+1]
                         + u_n[k+Nz] + u_n[k+1] - 6*u_n[k]);
            k += offset_r;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[k-offset_r]
                         + u_n[k-Nz-1] + u_n[k+Nz+1]
                         + u_n[k+Nz] + u_n[k+1] - 6*u_n[k]);
        }
        // z = Nz, r = 0 or Nx
        for (size_t c=1; c<Ny; c++) {
            k = c * (Ny+1) + Nz;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k+offset_r] + u_n[k+dim1]
                         + u_n[k-Nz-1] + u_n[k+Nz+1]
                         + u_n[k-1] + u_n[k-Nz] - 6*u_n[k]);
            k += offset_r;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                           u_n[k-dim1] + u_n[k-offset_r]
                         + u_n[k-Nz-1] + u_n[k+Nz+1]
                         + u_n[k-1] + u_n[k-Nz] - 6*u_n[k]);
        }
        // r=0, c=0, z=0
         u_np1[0] = 2 * u_n[0] - u_nm1[0] - gamma_t * (u_n[0]-u_nm1[0]) + alpha2[0] * (
                           u_n[offset_r] + u_n[dim1]
                         + u_n[offset_c] + u_n[Nz+1]
                         + u_n[Nz] + u_n[1] - 6*u_n[k]);
        // r=Nx, c=0, z=0
        k = offset_r;
        u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2 * (
                           u_n[k-dim1] + u_n[0]
                         + u_n[k+offset_c] + u_n[k+Nz+1]
                         + u_n[k+Nz] + u_n[k+1] - 6*u_n[k]);
        // r=0, c=Ny, z=0
        k = offset_c;
        u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2 * (
                           u_n[k+offset_r] + u_n[k+dim1]
                         + u_n[k-Nz-1] + u_n[0]
                         + u_n[k+Nz] + u_n[k+1] - 6*u_n[k]);
        // r=Nx, c=Ny, z=0
        k = offset_r + offset_c;
        u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2 * (
                           u_n[k-dim1] + u_n[offset_c]
                         + u_n[k-Nz-1] + u_n[offset_r]
                         + u_n[k+Nz] + u_n[k+1] - 6*u_n[k]);
        // r=0, c=0, z=Nz
        k = Nz;
        u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2 * (
                           u_n[k+offset_r] + u_n[k+dim1]
                         + u_n[k+offset_c] + u_n[k+Nz+1]
                         + u_n[k-1] + u_n[0] - 6*u_n[k]);
        // r=0, c=Ny, z=Nz
        k = offset_c+Nz;
        u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2 * (
                           u_n[k+offset_r] + u_n[k+dim1]
                         + u_n[k-Nz-1] + u_n[Nz]
                         + u_n[k-1] + u_n[offset_c] - 6*u_n[k]);
        // r=Nx, c=0, z=Nz
        k = offset_r + Nz;
        u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2 * (
                           u_n[k-dim1] + u_n[Nz]
                         + u_n[k+offset_c] + u_n[k+Nz+1]
                         + u_n[k-1] + u_n[offset_r] - 6*u_n[k]);
        // r=Nx, c=Ny, z=Nz
        k = offset_r + offset_c + Nz;
        u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2 * (
                           u_n[k-dim1] + u_n[k-offset_r]
                         + u_n[k-Nz-1] + u_n[k-offset_c]
                         + u_n[k-1] + u_n[k-Nz] - 6*u_n[k]);
        default:
            std::cout << "Invalid boundary condition option!" << std::endl;
    }
}

template<typename Real>
void WaveEquation<Real>::update_2d()
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
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
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
            offset_r = Nx * (Ny+1);
            Real kappa;
            Real delta = dt / dh;
            for (size_t c=0; c<Ny+1; c++) {
                k = offset_r + c;
                kappa    = (vp[c]*delta - 1) / (vp[c]*delta + 1);
                u_np1[c] = u_n[Ny+1+c] + kappa * (u_np1[Ny+1+c] - u_n[c]);
                kappa    = (vp[k]*delta - 1) / (vp[k]*delta + 1);
                u_np1[k] = u_n[k-Ny-1] + kappa * (u_np1[k-Ny-1] - u_n[k]);
            }
            for (size_t r=0; r<Nx+1; r++) {
                k = r * (Ny+1);
                kappa    = (vp[k]*delta - 1) / (vp[k]*delta + 1);
                u_np1[k] = u_n[k+1] + kappa * (u_np1[k+1] - u_n[k]);
                k = k + Ny;
                kappa    = (vp[k]*delta - 1) / (vp[k]*delta + 1);
                u_np1[k] = u_n[k-1] + kappa * (u_np1[k-1] - u_n[k]);
            }
            break;
        }
        /*case 3: // Nuemann bound cond
        {
            // r,c = 0,0
            u_np1[0]  = u_n[0] + init_cond*(u_n[0]-u_nm1[0])  + 2*alpha2[0]*(u_n[Ny+1] - 2*u_n[0] + u_n[1]);
            // r,c = 0,N_y
            u_np1[Ny] = u_n[Ny] + init_cond*(u_n[Ny]-u_nm1[Ny]) + 2*alpha2[Ny]*(u_n[2*Ny+1] - 2*u_n[Ny] + u_n[Ny-1]);
            // r,c = Nx,0
            k = Nx * (Ny+1);
            u_np1[k]  = u_n[k] + init_cond*(u_n[k]-u_nm1[k]) + 2*alpha2[k]*(u_n[k-Ny-1] - 2*u_n[k] + u_n[k+1]);
            // r,c = Nx,Ny
            k = offset_r + Ny;
            u_np1[k]  = u_n[k] + init_cond*(u_n[k]-u_nm1[k]) + 2*alpha2[k]*(u_n[k-1] - 2*u_n[k] + u_n[k-1]);
            // r = 0
            for (size_t c=1; c<Ny; c++) {
                u_np1[c] = u_n[c] + init_cond*(u_n[c]-u_nm1[c]) + alpha2[c]*(2*u_n[Ny+1+c] - 4*u_n[c] + u_n[c+1] + u_n[c-1]);
            }
            // c = 0
            for (size_t r=1; r<Nx; r++) {
                k = r * (Ny+1);
                u_np1[k] = u_n[k] + init_cond*(u_n[k]-u_nm1[k]) + alpha2[k]*(u_n[k+Ny+1] - 4*u_n[k] + u_n[k-Ny-1] + 2*u_n[k+1]);
            }
            // r = Nx
            offset_r = Nx * (Ny+1);
            for (size_t c=1; c<Ny; c++) {
                k = offset_r + c;
                u_np1[k] = u_n[k] + init_cond*(u_n[k]-u_nm1[k]) + alpha2[k]*(2*u_n[k-Ny-1] - 4*u_n[k] + u_n[k+1] + u_n[k-1]);
            }
            // c = Ny
            for (size_t r=1; r<Nx; r++) {
                offset_r = r * (Ny+1);
                k        = offset_r + Ny;
                u_np1[k] = u_n[k] + init_cond*(u_n[k]-u_nm1[k]) + alpha2[k]*(u_n[k+Ny+1] - 4*u_n[k] + u_n[k-Ny-1] + 2*u_n[k-1]);
            }
            break;
        }*/
        case 3: // Periodic cond
        {
            // r = 0
            offset_r = Nx * (Ny+1);
            for (size_t c=1; c<Ny; c++) {
                u_np1[c] = 2 * u_n[c] - u_nm1[c] - gamma_t * (u_n[c]-u_nm1[c]) + alpha2[c] * (
                       u_n[offset_r+c] + u_n[c+Ny+1] + u_n[c-1] + u_n[c+1] - 4*u_n[c]);
            }
            // c = 0
            for (size_t r=1; r<Nx; r++) {
                k = r * (Ny+1);
                u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                       u_n[k-Ny-1] + u_n[k+Ny+1] + u_n[k+Ny] + u_n[k+1] - 4*u_n[k]); 
            }
            // left corner
            u_np1[0] = 2 * u_n[0] - u_nm1[0] - gamma_t * (u_n[0]-u_nm1[0]) + alpha2[0] * (
                       u_n[offset_r] + u_n[Ny+1] + u_n[Ny] + u_n[1] - 4*u_n[0]); 
            k = offset_r;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                       u_n[k-Ny-1] + u_n[0] + u_n[k+Ny] + u_n[k+1] - 4*u_n[k]);
            // r = Nx 
            for (size_t c=1; c<Ny; c++) {
                k = c+offset_r;
                u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                       u_n[k-Ny-1] + u_n[c] + u_n[k-1] + u_n[k+1] - 4*u_n[k]); 
            }
            // c = Ny
            for (size_t r=1; r<Nx; r++) {
                k = r * (Ny+1) + Ny;
                u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                       u_n[k-Ny-1] + u_n[k+Ny+1] + u_n[k-1] + u_n[k-Ny] - 4*u_n[k]); 
            }
            // right corner
            u_np1[Ny] = 2 * u_n[Ny] - u_nm1[Ny] - gamma_t * (u_n[Ny]-u_nm1[Ny]) + alpha2[k] * (
                       u_n[0] + u_n[Ny-1] + u_n[2*Ny+1] + u_n[offset_r+Ny] - 4*u_n[Ny]);
            k = offset_r + Ny;
            u_np1[k] = 2 * u_n[k] - u_nm1[k] - gamma_t * (u_n[k]-u_nm1[k]) + alpha2[k] * (
                       u_n[k-Ny-1] + u_n[offset_r] + u_n[Ny] + u_n[k-1] - 4*u_n[k]); 
            break;
        }
        default:
            std::cout << "Invalid boundary condition option!" << std::endl;
    }
}
