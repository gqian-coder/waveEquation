#ifndef WAVEEQUATION_HPP
#define WAVEEQUATION_HPP 

template <typename Real> class WaveEquation {
public:
    //! Constructor.
    //!
    //!\param Nx and Ny simulation domain area
    //!\param dt and dh time and space grid spacings 
    //!\param C wave speed
    //!\param gamma dissipation rate
    //!\param use_condition Simulation boundary condiction
    //! use_condiction == 1: Dirichlet bound cond
    //! use_condiction == 2: Mur bound cond
    //! use_condiction == 3: Nuemann bound cond
    WaveEquation(const std::size_t Nx,
                 const std::size_t Ny,
                 const std::size_t Nz,
                 Real dt,
                 Real dh,
                 Real c,
                 Real gamma,
                 int use_condition);

    //! Update the simulation by one time tick.
    void update_2d(bool init_cond);
    void update_3d();

    // initial current ts data values
    void init_u_n(Real *data);

    // initial future ts data values
    void init_u_np1(Real *data);

    //! Return the size in bytes of the data
    std::size_t size() const;

    //! Vector array u_{i,j}^{n+1}
    std::vector<Real> u_np1;

private:
    //! Vector array u_{i,j}^{n-1} 
    std::vector<Real> u_nm1;
    //! Vector array u_{i,j}^{n}
    std::vector<Real> u_n; 
    std::size_t Nx;
    std::size_t Ny;
    std::size_t Nz;
    Real dt;
    Real dh;
    Real C; 
    Real gamma;
    Real alpha;
    Real alpha2;
    int use_condition;
};

#include "waveEquation.tpp"
#endif
