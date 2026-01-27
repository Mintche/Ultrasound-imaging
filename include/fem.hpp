#ifndef FEM_HPP
#define FEM_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include "math.hpp" // Ta classe matrice
#include "mesh.hpp" // La classe maillage 

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
    // Remplir le vecteur phi avec les valeurs des 6 fonctions de base au point (x,y)
    // Ordre des noeuds : 0,1,2 (sommets) puis 3,4,5 (milieux)
    static void evaluate_shape_functions(double x, double y, std::vector<double>& phi) {
        // Astuce : utilise les coordonnées barycentriques
        // L1 = 1 - x - y;
        // L2 = x;
        // L3 = y;

        // TODO: Implémenter les formules pour les sommets (phi[0] à phi[2])
        // ex: phi[0] = L1 * (2 * L1 - 1);

        // TODO: Implémenter les formules pour les milieux (phi[3] à phi[5])
        // ex: phi[3] = 4 * L1 * L2;
    }

    // -------------------------------------------------------------------------
    // 2. Gradients des fonctions de base
    // -------------------------------------------------------------------------
    // Remplir le vecteur grads avec les dérivées (d/dx, d/dy)
    static void evaluate_gradients(double x, double y, std::vector<std::pair<double, double>>& grads) {
        // TODO: Calculer les dérivées partielles par rapport à x et y
        // Attention : dL1/dx = -1, dL2/dx = 1, etc.
        
        // ex pour le noeud 0 :
        // grads[0].first  = ... (dérivée / x)
        // grads[0].second = ... (dérivée / y)
    }

    // -------------------------------------------------------------------------
    // 3. Quadrature de Gauss
    // -------------------------------------------------------------------------
    // Voir le TABLEAU page 5 du PDF
    static std::vector<QuadraturePoint> get_quadrature_points() {
        std::vector<QuadraturePoint> qp;
        
        // TODO: Remplir le vecteur avec les 7 points de Gauss.
        // Attention aux symétries indiquées dans le tableau (s1, s1), (s1, s3), etc.
        // Attention : Vérifier si les poids doivent être multipliés par 0.5 (aire du triangle ref)
        
        return qp;
    }

    // -------------------------------------------------------------------------
    // 4. Assemblage de la Matrice de Rigidité (Stiffness) A
    // -------------------------------------------------------------------------
    // A_ij = Integrale( grad(wi) . grad(wj) )
    static void assemble_stiffness(const Mesh& mesh, ProfileMatrix<complexe>& A) {
        // Récupérer les points de quadrature
        auto qp = get_quadrature_points();
        
        // Variables locales pour éviter les allocations dans la boucle
        std::vector<std::pair<double, double>> dphi_ref(6);
        
        // Boucle sur tous les triangles du maillage
        for (const auto& tri : mesh.triangles) {
            
            // TODO 1: Récupérer les coordonnées des 3 sommets du triangle (p0, p1, p2)
            
            // TODO 2: Calculer la Matrice Jacobienne J du passage (Reference -> Réel)
            // J = [ x1-x0  x2-x0 ]
            //     [ y1-y0  y2-y0 ]
            
            // TODO 3: Calculer le déterminant (detJ) et l'inverse de J
            // On a besoin de l'inverse pour transformer les gradients : 
            // Grad_reel = J^(-T) * Grad_ref
            
            // Boucle sur les points de quadrature
            for (const auto& q : qp) {
                // Evaluer les gradients sur le triangle de référence
                evaluate_gradients(q.x, q.y, dphi_ref);

                // Double boucle sur les fonctions de base (i et j de 0 à 5)
                for (int i = 0; i < 6; ++i) {
                    for (int j = 0; j < 6; ++j) {
                        
                        // TODO 4: Transformer les gradients ref en gradients réels
                        // dphi_dx_reel = (invJ_11 * dphi_dx_ref + invJ_12 * dphi_dy_ref)
                        
                        // TODO 5: Faire le produit scalaire grad_i . grad_j
                        
                        // TODO 6: Ajouter à la matrice globale A
                        // A(noeud_global_i, noeud_global_j) += val * poids * |detJ|
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // 5. Assemblage de la Matrice de Masse B
    // -------------------------------------------------------------------------
    // B_ij = Integrale( k^2 * wi * wj )
    static void assemble_mass(const Mesh& mesh, ProfileMatrix<complexe>& B, 
                              double k0, double k_d_val) {
        auto qp = get_quadrature_points();
        std::vector<double> phi(6);

        for (const auto& tri : mesh.triangles) {
            // TODO 1: Calculer le Jacobien (comme pour A) pour avoir detJ
            
            // TODO 2: Déterminer la valeur de k^2 pour ce triangle
            // Si tri.physical_tag indique un défaut -> utiliser k_d_val
            // Sinon -> utiliser k0
            // Voir section 1.2 "Modélisation des défauts"

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