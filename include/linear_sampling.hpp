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

class LinearSampling {
public:

    // -------------------------------------------------------------------------
    // Calcul de la composante du champs diffracté sur un bord 
    // -------------------------------------------------------------------------

    static void compute_boundary_u_s(const MeshP2& mesh, int n_mode,int tag_left,int tag_right, double direction, double L, double k0, double h, vector<complexe>& u_s, const vector<complexe>& u_n) {

        for (const auto& nodes : mesh.nodes) {
            if (nodes.ref == tag_left || nodes.ref == tag_right){
                double c_n_val = Fem::evaluate_c_1d(nodes.y,h,n_mode);
                complexe b_n = Fem::compute_beta(k0,h,n_mode);
                u_s[nodes.id] = u_n[nodes.id] - c_n_val * exp(direction * complexe(0., 1.) * b_n * nodes.x);
            }
        }
    }

    // -------------------------------------------------------------------------
    // Projection du champ diffracté n sur le mode m
    // -------------------------------------------------------------------------

    static void compute_projection(const vector<complexe>& u_s, const FullMatrix<complexe>& E_plus, const FullMatrix<complexe>& E_minus, vector<complexe>& U_proj_plus, vector<complexe>& U_proj_minus) {

        U_proj_plus = E_plus.transpose() * u_s; // Projection sur les modes de droite
        U_proj_minus = E_minus.transpose() * u_s; // Projection sur les modes de gauche

    }
    
    // -------------------------------------------------------------------------
    // Calcul de F
    // -------------------------------------------------------------------------

    static FullMatrix<complexe> compute_F(const FullMatrix<complexe>& S_LL, 
                                          const FullMatrix<complexe>& S_RL, 
                                          const FullMatrix<complexe>& S_LR, 
                                          const FullMatrix<complexe>& S_RR, 
                                          int N_modes, double k0, double h, double L) {

        FullMatrix<complexe> F(2*N_modes, 2*N_modes);

        for (int n = 0; n < N_modes; ++n) {
            
            // Calcul du terme de phase et de normalisation (Eq. 231)
            complexe beta_n = Fem::compute_beta(k0, h, n);
            complexe i_beta = complexe(0., 1.) * beta_n;
            
            // Constante = e^(i * beta * L) / (i * beta)
            complexe constante = exp(i_beta * L) / i_beta;

            for (int m = 0; m < N_modes; ++m) {
                // Remplissage selon Eq (7) et définitions (231-235)
                
                // Bloc F++ (Haut-Gauche) : Utilise (Un-)+ => Source Droite, Mesure Droite
                F(m, n) = constante * S_RR(n, m); 
                
                // Bloc F+- (Haut-Droite) : Utilise (Un+)+ => Source Gauche, Mesure Droite
                F(m, n + N_modes) = constante * S_RL(n, m);
                
                // Bloc F-+ (Bas-Gauche)  : Utilise (Un-)- => Source Droite, Mesure Gauche
                F(m + N_modes, n) = constante * S_LR(n, m);
                
                // Bloc F-- (Bas-Droite)  : Utilise (Un+)- => Source Gauche, Mesure Gauche
                F(m + N_modes, n + N_modes) = constante * S_LL(n, m);
            }
        }
        return F;
    }

    // -------------------------------------------------------------------------
    // Calcul de F++ (back scattering)
    // -------------------------------------------------------------------------

    static FullMatrix<complexe> compute_F_pp(const FullMatrix<complexe>& S_m_p, int N_modes, double k0, double h, double L) {

        FullMatrix<complexe> F_pp(2*N_modes,2*N_modes);

        for (int n = 0; n < N_modes; ++n) {

            complexe exposant = complexe(0.,1.) * Fem::compute_beta(k0, h, n) * L; // exp(i*beta_n*L)
            complexe constante = exp(exposant) / (complexe(0.,1.) * Fem::compute_beta(k0, h, n)); // exp(i*beta_n*L) / (i*beta_n)

            for (int m = 0; m < N_modes; ++m) {
                F_pp(m, n) = constante * S_m_p(n,m); //F++ U-+
            }
        }
        return F_pp;
    }

    // -------------------------------------------------------------------------
    // Calcul de G 
    // -------------------------------------------------------------------------

    static vector<complexe> assemble_Gz(const MeshP2& mesh, int n_modes, double z1, double z2, double x_min, double x_max, double k0, double h){
        
        vector<complexe> Gz(2*n_modes);
        
        for (int m = 0; m < n_modes;m++){
            
            complexe denominateur = complexe(0.,1.) * Fem::compute_beta(k0, h, m);
            complexe c_val = Fem::evaluate_c_1d(z2, h, m); // Attention: suppose y dans [0, h]
            
            // Phase relative au bord gauche (x_min) et droit (x_max)
            Gz[m] = exp(complexe(0.,1.) * Fem::compute_beta(k0, h, m) * (x_max - z1)) * c_val / denominateur;
            Gz[m + n_modes] = exp(complexe(0.,1.) * Fem::compute_beta(k0, h, m) * (z1 - x_min)) * c_val / denominateur;
        }
        return Gz;
    }

    // -------------------------------------------------------------------------
    // Bruit Gaussien
    // -------------------------------------------------------------------------

    static void add_gaussian_noise(vector<complexe>& data, double noise_level) {
        std::random_device rd;
        std::mt19937 generator(rd());
        std::normal_distribution<double> distribution(0.0, noise_level);

        for (auto& val : data) {
            double noise_real = distribution(generator);
            double noise_imag = distribution(generator);
            val += complexe(noise_real, noise_imag);
        }
    }

    
};

#endif // LINEAR_SAMPLING_HPP