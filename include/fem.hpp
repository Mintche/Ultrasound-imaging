#ifndef FEM_HPP
#define FEM_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include "math.hpp" // Ta classe matrice
#include "mesh.hpp" // La classe maillage

using namespace std;
using namespace usim;

// Structure pour stocker un point d'intégration (Gauss)
struct QuadraturePoint {
    double x, y; // Coordonnées sur le triangle de référence
    double w;    // Poids associé
};

class Fem {
public:
    // -------------------------------------------------------------------------
    // 1. Fonctions de base P2 (Shape Functions) sur le triangle de référence
    // -------------------------------------------------------------------------

    // Ordre des noeuds : 0,1,2 (sommets) puis 3,4,5 (milieux)
    static void evaluate_shape_functions(double x, double y, vector<double>& phi) {

        double lambda1 = 1.0 - x - y;
        double lambda2 = x;
        double lambda3 = y;

        phi[0] = lambda1 * (2.0 * lambda1 - 1.0); // Sommet 0
        phi[1] = lambda2 * (2.0 * lambda2 - 1.0); // Sommet 1
        phi[2] = lambda3 * (2.0 * lambda3 - 1.0); // Sommet 2

        phi[3] = 4.0 * lambda1 * lambda2; // Milieu arête 0-1
        phi[4] = 4.0 * lambda2 * lambda3; // Milieu arête 1-2
        phi[5] = 4.0 * lambda1 * lambda3; // Milieu arête 0-2
    }

    // -------------------------------------------------------------------------
    // 2. Gradients des fonctions de base
    // -------------------------------------------------------------------------

    static void evaluate_gradients(double x, double y, FullMatrix<double>& grads) {

        grads(0,0) = 4.0*x + 4.0*y - 3.0;
        grads(0,1) = 4.0*x + 4.0*y - 3.0;

        grads(1,0) = 4.0*x - 1.0;
        grads(1,1) = 0.0;

        grads(2,0) = 0.0;
        grads(2,1) = 4.0*y - 1.0;

        grads(3,0) = 4.0 - 8.0*x - 4.0*y;
        grads(3,1) = -4.0*x;

        grads(4,0) = 4.0*y;
        grads(4,1) = 4.0*x;

        grads(5,0) = -4.0*y;
        grads(5,1) = 4.0 - 4.0*x - 8.0*y;
    }

    // -------------------------------------------------------------------------
    // 3. Quadrature de Gauss
    // -------------------------------------------------------------------------

    static vector<QuadraturePoint> get_quadrature_points() {

        double s0 = 1/3.;
        double s1 = (6-sqrt(15))/21.;
        double s2 = (6+sqrt(15))/21.;
        double s3 = (9+2*sqrt(15))/21.;
        double s4 = (9-2*sqrt(15))/21.;

        double eta0 = 9/80.;
        double eta1 = (155-sqrt(15))/2400.;
        double eta2 = (155+sqrt(15))/2400.;

        vector<QuadraturePoint> qp(7);

        qp[0].x = s0;
        qp[0].y = s0;
        qp[0].w = eta0;

        qp[1].x = s1;
        qp[1].y = s1;
        qp[1].w = eta1;

        qp[2].x = s1;
        qp[2].y = s3;
        qp[2].w = eta1;

        qp[3].x = s3;
        qp[3].y = s1;
        qp[3].w = eta1;

        qp[4].x = s2;
        qp[4].y = s2;
        qp[4].w = eta2;

        qp[5].x = s2;
        qp[5].y = s4;
        qp[5].w = eta2;

        qp[6].x = s4;
        qp[6].y = s2;
        qp[6].w = eta2;
        
        return qp;

    }

    // -------------------------------------------------------------------------
    // 4. Calcul du profil de la matrice
    // -------------------------------------------------------------------------

    static vector<size_t> compute_profile(const MeshP2& mesh) {
        size_t ndof = mesh.ndof();
        vector<size_t> p(ndof);
        
        // Initialisation : p[i] = i (au moins la diagonale)
        for(size_t i=0; i<ndof; ++i) p[i] = i;

        // Parcours des éléments pour mettre à jour la largeur de bande
        for(const auto& tri : mesh.triangles) {
            for(int i=0; i<6; ++i) {
                int u = tri.node_ids[i]; // Ligne potentielle
                for(int j=0; j<6; ++j) {
                    int v = tri.node_ids[j]; // Colonne potentielle
                    // Si v < p[u], on élargit le profil
                    if (v < static_cast<int>(p[u])) {
                        p[u] = v;
                    }
                }
            }
        }
        return p;
    }

    // -------------------------------------------------------------------------
    // 5. Assemblage de la Matrice de Rigidité (Stiffness) A
    // -------------------------------------------------------------------------
    // A_ij = Integrale( grad(wi) . grad(wj) )
    static void assemble_stiffness(const MeshP2& mesh, ProfileMatrix<complexe>& A) {
        // Récupérer les points de quadrature
        vector<QuadraturePoint> qp = get_quadrature_points();
        
        // Variables locales pour éviter les allocations dans la boucle
        FullMatrix<double> dphi_ref(6,2);
        
        // Boucle sur tous les triangles du maillage
        for (const auto& tri : mesh.triangles) {
            
            // Coordonnées des 3 sommets du triangle (p0, p1, p2)

            Point2D p0 = mesh.nodes[tri.node_ids[0]];
            Point2D p1 = mesh.nodes[tri.node_ids[1]];
            Point2D p2 = mesh.nodes[tri.node_ids[2]];

            //Passage (Reference -> Réel) F(S) = Bl*S + bl avec bl = S0, Bl = [S1-S0,S2-S0]

            FullMatrix<double> Jac(2,2);

            Jac(0,0) = p1.x-p0.x;
            Jac(0,1) = p2.x-p0.x;

            Jac(1,0) = p1.y-p0.y;
            Jac(1,1) = p2.y-p0.y;

            //detJ et l'inverse de J

            double detJac = (p1.x-p0.x)*(p2.y-p0.y)-(p1.y-p0.y)*(p2.x-p0.x);

            FullMatrix<double> invJac = Jac.inverse();

            // Matrice de rigidité élémentaire (6x6) pour ce triangle
            FullMatrix<double> A_elem(6, 6); // Initialisée à 0
            
            // Boucle sur les points de quadrature
            for (const auto& q : qp) {
                // Evaluer les gradients sur le triangle de référence
                evaluate_gradients(q.x, q.y, dphi_ref);

                // 1. Gradients réels : G_real = dphi_ref * invJac
                // dphi_ref (6x2) * invJac (2x2) -> G_real (6x2)
                FullMatrix<double> G_real = dphi_ref * invJac;

                // 2. Contribution locale : G_real * G_real^T
                // Cela calcule tous les produits scalaires grad_i . grad_j d'un coup
                FullMatrix<double> contrib = G_real * G_real.transpose();

                // 3. Accumulation pondérée dans K_elem
                double weight = q.w * abs(detJac);
                
                // On ajoute manuellement car pas d'opérateur K_elem += contrib * scalar
                for(int i=0; i<6; ++i) {
                    for(int j=0; j<6; ++j) {
                        A_elem(i, j) += contrib(i, j) * weight;
                    }
                }
            }
            
            // 4. Assemblage dans la matrice
            for(int i=0; i<6; ++i) {
                int I = tri.node_ids[i];
                for(int j=0; j<6; ++j) {
                    int J = tri.node_ids[j];
                    A(I, J) += complexe(A_elem(i, j), 0.0);
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // 6. Assemblage de la Matrice de Masse B
    // -------------------------------------------------------------------------

    // B_ij = Integrale( k^2 * wi * wj )
    static void assemble_mass(const MeshP2& mesh, ProfileMatrix<complexe>& B, 
                              double k0, double k_d_val) {
        auto qp = get_quadrature_points();
        std::vector<double> phi(6);

        for (const auto& tri : mesh.triangles) {

            // Coordonnées des 3 sommets du triangle (p0, p1, p2)

            Point2D p0 = mesh.nodes[tri.node_ids[0]];
            Point2D p1 = mesh.nodes[tri.node_ids[1]];
            Point2D p2 = mesh.nodes[tri.node_ids[2]];

            //Passage (Reference -> Réel) F(S) = Bl*S + bl avec bl = S0, Bl = [S1-S0,S2-S0]

            FullMatrix<double> Jac(2,2);

            Jac(0,0) = p1.x-p0.x;
            Jac(0,1) = p2.x-p0.x;

            Jac(1,0) = p1.y-p0.y;
            Jac(1,1) = p2.y-p0.y;

            //detJ et l'inverse de J

            double detJac = (p1.x-p0.x)*(p2.y-p0.y)-(p1.y-p0.y)*(p2.x-p0.x);

            FullMatrix<double> invJac = Jac.inverse();
            
            double k = (tri.ref == 0 ) ? k0 : k_d_val;

            // Voir section 1.2 "Modélisation des défauts"

            // Boucle quadrature
            for (const auto& q : qp) {
                evaluate_shape_functions(q.x, q.y, phi);

                double weight = q.w * abs(detJac);

                // Double boucle sur les fonctions de base
                for (int i = 0; i < 6; ++i) {
                    for (int j = 0; j < 6; ++j) {
                        
                        double val = phi[i] * phi[j] * weight * (k * k);
                        
                        B(tri.node_ids[i],tri.node_ids[j]) += val;
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // 7. Quadrature 1D (Gauss-Legendre) pour les bords
    // -------------------------------------------------------------------------
    static vector<QuadraturePoint> get_quadrature_points_1d() {
        vector<QuadraturePoint> qp(3);
        double sqrt35 = sqrt(3/5.); // sqrt(3/5)
        
        qp[0].x = -sqrt35; qp[0].w = 5.0/9.0;
        qp[1].x = 0.0;     qp[1].w = 8.0/9.0;
        qp[2].x = sqrt35;  qp[2].w = 5.0/9.0;
        
        // .y n'est pas utilisé en 1D
        return qp;
    }

    // -------------------------------------------------------------------------
    // 8. Fonctions de forme 1D P2 sur [-1, 1] et calcul de c_n
    // -------------------------------------------------------------------------
    // phi[0] : t=-1 (gauche), phi[1] : t=1 (droite), phi[2] : t=0 (milieu)
    static void evaluate_shape_functions_1d(double t, vector<double>& phi) {
        phi[0] = 0.5 * t * (t - 1.0); // Sommet gauche
        phi[1] = 0.5 * t * (t + 1.0); // Sommet droite
        phi[2] = 1.0 - t * t;         // Milieu
    }

    static double evaluate_c_1d(double x, double h, int n) {
        return (2.0 / h)*cos(n*M_PI * x / h);        
    }

    // -------------------------------------------------------------------------
    // 9. Assemblage générique d'une matrice de masse surfacique (Bord)
    // -------------------------------------------------------------------------
    // Calcule Integrale_Gamma ( coeff * u * v ) sur les arêtes ayant le tag 'boundary_tag'
    static void assemble_E(const MeshP2& mesh, FullMatrix<complex<double>>& E, int N, int boundary_tag) {
        auto qp = get_quadrature_points_1d();
        vector<double> phi_1d(3);
        double h = mesh.Ly; // Hauteur du domaine 
        
        for (const auto& tri : mesh.triangles) {
            for (int edge_i = 0; edge_i < 3; ++edge_i) {
                // Si l'arête courante possède le tag recherché
                if (tri.edge_ref[edge_i] == boundary_tag) {
                    // Calcul de F' pour cette arête
                    int idx_A = edge_i;
                    int idx_B = (edge_i + 1) % 3;
                    int idx_M = edge_i + 3;

                    int nodes_global[3] = {
                        tri.node_ids[idx_A],
                        tri.node_ids[idx_B], 
                        tri.node_ids[idx_M]
                    };  
                    Point2D A = mesh.nodes[nodes_global[0]];
                    Point2D B = mesh.nodes[nodes_global[1]];
                    double edge_length = sqrt((B.x - A.x)*(B.x - A.x) + (B.y - A.y)*(B.y - A.y));
                    double detJac = edge_length / 2.0;

                    for (const auto& q : qp) {
                        double t = q.x; // Coordonnée sur l'arête de référence [-1,1]
                        double weight = q.w;
                        evaluate_shape_functions_1d(t, phi_1d);

                        for (int i = 0; i < 3 ; i++) {
                            for (int n = 0; n < N; n++) {
                                double c_n = evaluate_c_1d(0.5 * ( (B.x + A.x) + t * (B.x - A.x) ), h, n);
                                complex<double> val = phi_1d[i] * phi_1d[n] * c_n * weight * detJac;
                                E(nodes_global[i], n) += val;           
                            }
                        }
                    } 
                }
            }
        }
    }

    static void assemble_D(const MeshP2& mesh, FullMatrix<complex<double>>& D, int boundary_tag, int N, int total_dofs, double k0) {
        int mini = min(total_dofs, N);
            for (int j =0; j<mini; j++) {
                D(j,j) = complex<double>(0.0,1.0)* compute_beta(k0, mesh.Ly, j);
            }       
    }

    static void assemble_T(FullMatrix<complex<double>>& T, FullMatrix<complex<double>>& D, FullMatrix<complex<double>>& E, int N, int boundary_tag){
        T = E * D * E.transpose();
    }
        



// -------------------------------------------------------------------------
// 10. Calcul de Beta_n
// -------------------------------------------------------------------------

    complexe compute_beta(double k0, double h, int n) {
        return sqrt(complexe(k0*k0 - (M_PI*n/h)*(M_PI*n/h), 0.0));
    }

};

#endif