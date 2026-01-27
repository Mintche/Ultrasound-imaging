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

        phi[0] = (1-x-y)*(1-x-2*y);
        phi[1] = x*2*(x-0.5);
        phi[2] = y*2*(y-0.5);

        phi[3] = 4*(1-x-y)*y;
        phi[4] = 4*(1-x-y)*x;
        phi[5] = 4*x*y;
    }

    // -------------------------------------------------------------------------
    // 2. Gradients des fonctions de base
    // -------------------------------------------------------------------------

    static void evaluate_gradients(double x, double y, FullMatrix<double>& grads) {

        grads(0,0) = 2*x+3*y-2;  // dérivée par rapport à x
        grads(0,1) = 2*x+4*y-3; // dérivée par rapport à y

        grads(1,0) = 4*x-1;
        grads(1,1) = 0;

        grads(2,0) = 4*y-1;
        grads(2,1) = 0;

        grads(3,0) = -4*y*x;
        grads(3,1) = 4*(1-x-2*y);

        grads(4,0) = 4*(1-y-2*x);
        grads(4,1) = -4*x*y;

        grads(5,0) = 4*y;
        grads(5,1) = 4*x;
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
                double weight = q.w * std::abs(detJac);
                
                // On ajoute manuellement car pas d'opérateur K_elem += contrib * scalar
                for(int i=0; i<6; ++i) {
                    for(int j=0; j<6; ++j) {
                        A_elem(i, j) += contrib(i, j) * weight;
                    }
                }
            }
            
            // 4. Assemblage dans la matrice globale A (une seule fois par triangle)
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
            // TODO 1: Calculer le Jacobien (comme pour A) pour avoir detJ
            
            // TODO 2: Déterminer la valeur de k^2 pour ce triangle
            // Si tri.physical_tag indique un défaut -> utiliser k_d_val
            // Sinon -> utiliser k0
            // Voir section 1.2 "Modélisation des défauts" [cite: 54]

            // Boucle quadrature
            for (const auto& q : qp) {
                evaluate_shape_functions(q.x, q.y, phi);

                // Double boucle sur les fonctions de base
                for (int i = 0; i < 6; ++i) {
                    for (int j = 0; j < 6; ++j) {
                        
                        // TODO 3: Calculer l'intégrale locale
                        // val = k^2 * phi[i] * phi[j]
                        
                        // TODO 4: Ajouter à la matrice B (avec le poids et detJ)
                    }
                }
            }
        }
    }
};

#endif