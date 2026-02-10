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
    static void compute_u_s(const MeshP2& mesh, int n_mode, int boundary_tag, double direction, double L,double k0, double h,vector<complexe>& u_s, vector<complexe>& u_n) {

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
                            complexe exposant = exp(direction * complexe(0.,1.) *Fem::compute_beta(k0, h, n_mode)*L); // exp(i*direction*dot((x,y), direction_vector))
                            complexe phi_n = Fem::evaluate_c_1d(y, h, n_mode)*exp(exposant); // Calcul de phi_n à ce point
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
    
    // -------------------------------------------------------------------------
    // Calcul de F
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Calcul de F++ (back scattering)
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Calcul de G 
    // -------------------------------------------------------------------------

    static FullMatrix<complexe> assemble_F_G(const MeshP2& mesh, int n_mode, double z1, double z2, double L,double k0, double h){
        FullMatrix<complexe> F_G(n_mode,2);
        for (int m = 0; m < n_mode;m++){
            F_G(m,1) = exp(complexe(0.,1.) *Fem::compute_beta(k0, h, n_mode)*(L+z1))*Fem::evaluate_c_1d(z2,h,m)/(complexe(0.,1.) *Fem::compute_beta(k0, h, n_mode)*(L+z1));
            F_G(m,2) = exp(complexe(0.,1.) *Fem::compute_beta(k0, h, n_mode)*(L-z1))*Fem::evaluate_c_1d(z2,h,m)/(complexe(0.,1.) *Fem::compute_beta(k0, h, n_mode)*(L+z1));
        }
        return F_G;
    }

    // -------------------------------------------------------------------------
    // Generation d'image matlab de log(1/||h||)
    // -------------------------------------------------------------------------
};
#endif // LINEAR_SAMPLING_HPP
