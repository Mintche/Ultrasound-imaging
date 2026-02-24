#ifndef LINEAR_SAMPLING_HPP
#define LINEAR_SAMPLING_HPP

#include <cmath>
#include <complex>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include <random>

#include "fem.hpp"

namespace LinearSampling {
    // -------------------------------------------------------------------------
    // Calcul de la composante du champs diffracté sur un bord 
    // -------------------------------------------------------------------------

    void compute_boundary_u_s(const MeshP2& mesh, int n_mode,int tag_left,int tag_right,
                            double direction, double L, double k0, double h, std::vector<complexe>& u_s, 
                            const std::vector<complexe>& u_n);

    // -------------------------------------------------------------------------
    // Projection du champ diffracté n sur le mode m
    // -------------------------------------------------------------------------

    void compute_projection(const std::vector<complexe>& u_s, const FullMatrix<complexe>& E_plus, 
                            const FullMatrix<complexe>& E_minus, std::vector<complexe>& U_proj_plus, 
                            std::vector<complexe>& U_proj_minus);
    
    // -------------------------------------------------------------------------
    // Calcul de F
    // -------------------------------------------------------------------------

    FullMatrix<complexe> compute_F(const FullMatrix<complexe>& S_LL, 
                                    const FullMatrix<complexe>& S_RL, 
                                    const FullMatrix<complexe>& S_LR, 
                                    const FullMatrix<complexe>& S_RR, 
                                    int N_modes, double k0, double h, double L);

    // -------------------------------------------------------------------------
    // Calcul de F++ (back scattering)
    // -------------------------------------------------------------------------

    FullMatrix<complexe> compute_F_pp(const FullMatrix<complexe>& S_m_p, int N_modes, double k0, double h, double L);

    // -------------------------------------------------------------------------
    // Calcul de G 
    // -------------------------------------------------------------------------

    std::vector<complexe> assemble_Gz(const MeshP2& mesh, int n_modes, double z1, double z2, double x_min, double x_max, double k0, double h);

    // -------------------------------------------------------------------------
    // Bruit Gaussien
    // -------------------------------------------------------------------------

    void add_gaussian_noise(std::vector<complexe>& data, double noise_level) ;
};

#endif // LINEAR_SAMPLING_HPP