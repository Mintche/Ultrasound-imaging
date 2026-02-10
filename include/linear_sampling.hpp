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

#include "fem.hpp"

class LinearSampling {
public:

    // -------------------------------------------------------------------------
    // Calcul de la composante du champs diffracté sur un bord 
    // -------------------------------------------------------------------------

    static void compute_u_s(const MeshP2& mesh, int n_modes, int boundary_tag, double direction, double L, double k0, double h, vector<complexe>& u_s, const vector<complexe>& u_n) {

        vector<int> visited(mesh.ndof(), 0);
        for (const auto& tri : mesh.triangles) {
            for (int edge_i = 0; edge_i < 3; ++edge_i) {
                if (tri.edge_ref[edge_i] == boundary_tag){
                    int idx_A = edge_i;
                    int idx_B = (edge_i + 1) % 3;
                    int idx_M = edge_i + 3;

                    // Indices globaux dans la matrice
                    int nodes_global[3] = {
                        tri.node_ids[idx_A],
                        tri.node_ids[idx_B], 
                        tri.node_ids[idx_M]
                    };

                    for (const auto& idx : nodes_global) {
                        if (!visited[idx]) {
                            visited[idx] = 1;
                            double x = mesh.nodes[idx].x;
                            double y = mesh.nodes[idx].y;

                            complexe exposant = direction * complexe(0.,1.) * Fem::compute_beta(k0, h, n_modes) * x; //on prend x pour le bord même si x = L
                            complexe phi_n = Fem::evaluate_c_1d(y, h, n_modes) * exp(exposant); 
                            
                            u_s[idx] = u_n[idx] - phi_n; // Contribution au champ diffracté sur le bord
                        }        
                    }        
                }   
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

    static FullMatrix<complexe> compute_F(const FullMatrix<complexe>& S_p_m, const FullMatrix<complexe>& S_m_p, const FullMatrix<complexe>& S_p_p, const FullMatrix<complexe>& S_m_m, int N_modes, double k0, double h, double L) {

        FullMatrix<complexe> F(2*N_modes,2*N_modes);

        for (int n = 0; n < N_modes; ++n) {

            complexe exposant = complexe(0.,1.) * Fem::compute_beta(k0, h, n) * L; // exp(i*beta_n*L)
            complexe constante = exp(exposant) / (complexe(0.,1.) * Fem::compute_beta(k0, h, n)); // exp(i*beta_n*L) / (i*beta_n)

            for (int m = 0; m < N_modes; ++m) {
                F(m, n) = constante * S_m_p(n,m);                               // F++ U-+ 
                F(m, n + N_modes) = constante * S_p_p(n,m);                     // F+- U++
                F(m + N_modes, n) = constante * S_m_m(n,m);                     // F-+ U--
                F(m + N_modes, n + N_modes) = constante * S_p_m(n,m);           // F-- U+-
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

    static FullMatrix<complexe> assemble_F_G(const MeshP2& mesh, int n_modes, double z1, double z2, double L,double k0, double h){
        
        FullMatrix<complexe> F_G(n_modes,2);
        
        for (int m = 0; m < n_modes;m++){
            
            complexe denominateur = complexe(0.,1.) * Fem::compute_beta(k0, h, m);
            
            F_G(m,0) = exp(complexe(0.,1.) * Fem::compute_beta(k0, h, m) * (L+z1)) * Fem::evaluate_c_1d(z2,h,m) / denominateur;
            F_G(m,1) = exp(complexe(0.,1.) * Fem::compute_beta(k0, h, m) * (L-z1)) * Fem::evaluate_c_1d(z2,h,m) / denominateur;
        }
        return F_G;
    }

    // -------------------------------------------------------------------------
    // Generation d'image matlab de log(1/||h||)
    // -------------------------------------------------------------------------
    
};

#endif // LINEAR_SAMPLING_HPP