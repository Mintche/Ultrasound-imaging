#ifndef FEM_HPP
#define FEM_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "math.hpp"
#include "mesh.hpp"

using namespace std;
using namespace usim;

// Structure pour stocker un point d'intégration (Gauss)
struct QuadraturePoint {
    double x, y; // Coordonnées sur le triangle de référence
    double w;    // Poids associé
};

namespace Fem {
    // -------------------------------------------------------------------------
    // Fonctions de base P2 (Shape Functions) sur le triangle de référence
    // -------------------------------------------------------------------------

    // Ordre des noeuds : 0,1,2 (sommets) puis 3,4,5 (milieux)
    void evaluate_shape_functions(double x, double y, vector<double>& phi) {

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
    // Gradients des fonctions de base
    // -------------------------------------------------------------------------

    void evaluate_gradients(double x, double y, FullMatrix<double>& grads) {

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
    // Quadrature de Gauss
    // -------------------------------------------------------------------------

    vector<QuadraturePoint> get_quadrature_points() {

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
    // Calcul du profil de la matrice
    // -------------------------------------------------------------------------

    vector<size_t> compute_profile(const MeshP2& mesh){
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

    vector<size_t> compute_profile_enhanced(const MeshP2& mesh, const vector<int>& boundary_tags){
        // 1. Profil standard basé sur les triangles
        vector<size_t> p = compute_profile(mesh);

        // 2. Élargissement du profil pour les frontières

        for (int tag : boundary_tags) {
            vector<int> boundary_nodes;
            // On parcourt les arêtes pour trouver les noeuds du bord 'tag'
            for(const auto& tri : mesh.triangles) {
                for(int i=0; i<3; ++i) {
                    if(tri.edge_ref[i] == tag) {
                        // Les 3 noeuds de l'arête (2 sommets + 1 milieu)
                        boundary_nodes.push_back(tri.node_ids[i]);
                        boundary_nodes.push_back(tri.node_ids[(i+1)%3]);
                        boundary_nodes.push_back(tri.node_ids[i+3]);
                    }
                }
            }
            // Suppression des doublons
            sort(boundary_nodes.begin(), boundary_nodes.end());
            boundary_nodes.erase(unique(boundary_nodes.begin(), boundary_nodes.end()), boundary_nodes.end());

            // Pour chaque paire de noeuds (u, v) sur ce bord, on met à jour le profil
            // car la matrice T va créer un coefficient non nul entre eux.
            for (int u : boundary_nodes) {
                for (int v : boundary_nodes) {
                    if (u > v) { // On ne stocke que le triangle inférieur
                        if (v < static_cast<int>(p[u])) {
                            p[u] = v;
                        }
                    }
                }
            }
        }
        return p;
    }

    // -------------------------------------------------------------------------
    // Assemblage de la Matrice de Rigidité A
    // -------------------------------------------------------------------------

    void A_matrix(const MeshP2& mesh, ProfileMatrix<complexe>& A, double factor = 1.0){
        vector<QuadraturePoint> qp = get_quadrature_points();

        FullMatrix<double> dphi_ref(6,2);
        for (const auto& tri : mesh.triangles) {

            Point2D p0 = mesh.nodes[tri.node_ids[0]];
            Point2D p1 = mesh.nodes[tri.node_ids[1]];
            Point2D p2 = mesh.nodes[tri.node_ids[2]];

            //Passage (Reference -> Réel) F(S) = Bl*S + bl avec bl = S0, Bl = [S1-S0,S2-S0]

            double J00 = p1.x-p0.x; double J01 = p2.x-p0.x;
            double J10 = p1.y-p0.y; double J11 = p2.y-p0.y;
            double detJac = J00*J11 - J01*J10;
            double invDet = 1.0 / detJac;
            
            double iJ00 =  J11 * invDet;
            double iJ01 = -J01 * invDet;
            double iJ10 = -J10 * invDet;
            double iJ11 =  J00 * invDet;

            // Boucle sur les points de quadrature
            for (const auto& q : qp) {
                evaluate_gradients(q.x, q.y, dphi_ref);
                double weight = q.w * abs(detJac) * factor;

                double G[6][2];
                for(int i=0; i<6; ++i) {
                    G[i][0] = dphi_ref(i,0)*iJ00 + dphi_ref(i,1)*iJ10;
                    G[i][1] = dphi_ref(i,0)*iJ01 + dphi_ref(i,1)*iJ11;
                }
                
                for(int i=0; i<6; ++i) {
                    for(int j=0; j<=i; ++j) {
                        double dot = G[i][0]*G[j][0] + G[i][1]*G[j][1];
                        A(tri.node_ids[i], tri.node_ids[j]) += complexe(dot * weight, 0.0);
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Assemblage de la Matrice de Masse B
    // -------------------------------------------------------------------------

    void B_matrix(const MeshP2& mesh, ProfileMatrix<complexe>& B, 
                              double k0, double k_d_val, double factor = -1.0){
        auto qp = get_quadrature_points();
        vector<double> phi(6);

        for (const auto& tri : mesh.triangles) {

            // Coordonnées des 3 sommets du triangle (p0, p1, p2)

            Point2D p0 = mesh.nodes[tri.node_ids[0]];
            Point2D p1 = mesh.nodes[tri.node_ids[1]];
            Point2D p2 = mesh.nodes[tri.node_ids[2]];

            // detJ (Jacobien géométrique)
            double detJac = (p1.x-p0.x)*(p2.y-p0.y)-(p1.y-p0.y)*(p2.x-p0.x);

            // CORRECTION : Utilisation de is_defect pour choisir le bon k
            double k = (tri.is_defect) ? k_d_val : k0;

            // Voir section 1.2 "Modélisation des défauts"

            // Boucle quadrature
            for (const auto& q : qp) {
                evaluate_shape_functions(q.x, q.y, phi);

                double weight = q.w * abs(detJac);

                // Double boucle sur les fonctions de base
                for (int i = 0; i < 6; ++i) {
                    for (int j = 0; j <= i; ++j) {
                        
                        double val = phi[i] * phi[j] * weight * (k * k);
                        
                        B(tri.node_ids[i],tri.node_ids[j]) += val * factor;
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Quadrature 1D (Gauss-Legendre) pour les bords
    // -------------------------------------------------------------------------
    vector<QuadraturePoint> get_quadrature_points_1d(){
        vector<QuadraturePoint> qp(3);
        double sqrt35 = sqrt(3/5.);
        
        qp[0].x = -sqrt35; qp[0].w = 5.0/9.0;
        qp[1].x = 0.0;     qp[1].w = 8.0/9.0;
        qp[2].x = sqrt35;  qp[2].w = 5.0/9.0;
        
        return qp;
    }

    // -------------------------------------------------------------------------
    // Fonctions de forme 1D P2 sur [-1, 1] et calcul de c_n
    // -------------------------------------------------------------------------
    // phi[0] : t=-1 (gauche), phi[1] : t=1 (droite), phi[2] : t=0 (milieu)
    void evaluate_shape_functions_1d(double t, vector<double>& phi){
        phi[0] = 0.5 * t * (t - 1.0); // Sommet gauche
        phi[1] = 0.5 * t * (t + 1.0); // Sommet droite
        phi[2] = 1.0 - t * t;         // Milieu
    }

    double evaluate_c_1d(double y, double h, int n){
        if (n == 0) return sqrt(1.0 / h);
        return sqrt(2.0 / h) * cos(n * M_PI * y / h);
    }

    // -------------------------------------------------------------------------
    // Calcul de Beta_n
    // -------------------------------------------------------------------------

    complexe compute_beta(double k0, double h, int n){
        return sqrt(complexe(k0*k0 - (M_PI*n/h)*(M_PI*n/h), 0.0));
    }


    // -------------------------------------------------------------------------
    // Assemblage matrice E
    // -------------------------------------------------------------------------

    FullMatrix<complexe> compute_E(const MeshP2& mesh, int N_modes, int boundary_tag, double k0) {
        
        int ndof = mesh.ndof();
        // E est une matrice (Ndof x N_modes)
        FullMatrix<complexe> E(ndof, N_modes); 

        auto qp = get_quadrature_points_1d(); // Quadrature de Gauss 1D
        vector<double> phi_1d(3);
        double h = mesh.Ly; // Hauteur du guide

        // Boucle sur tous les triangles pour trouver les arêtes du bord
        for (const auto& tri : mesh.triangles) {
            for (int edge_i = 0; edge_i < 3; ++edge_i) {
                
                // Si l'arête est sur la frontière demandée (ex: tag 1 pour gauche, 2 pour droite)
                if (tri.edge_ref[edge_i] == boundary_tag) {
                    
                    // Indices locaux des noeuds de l'arête (A, B) et milieu (M)
                    int idx_A = edge_i;
                    int idx_B = (edge_i + 1) % 3;
                    int idx_M = edge_i + 3;

                    // Indices globaux dans la matrice
                    int nodes_global[3] = {
                        tri.node_ids[idx_A],
                        tri.node_ids[idx_B], 
                        tri.node_ids[idx_M]
                    };

                    // Coordonnées physiques pour calculer la longueur (Jacobien)
                    Point2D A = mesh.nodes[nodes_global[0]];
                    Point2D B = mesh.nodes[nodes_global[1]];
                    
                    // Longueur de l'arête (segment vertical)
                    double edge_length = sqrt(pow(B.x - A.x, 2) + pow(B.y - A.y, 2));
                    double detJac = edge_length / 2.0;

                    // Boucle sur les points d'intégration
                    for (const auto& q : qp) {
                        double t = q.x; // Coordonnée ref [-1, 1]
                        double w = q.w;
                        
                        // 1. Fonctions de forme 1D au point t
                        evaluate_shape_functions_1d(t, phi_1d);

                        // 2. Coordonnée physique Y au point d'intégration

                        double y_phys = 0.5 * ((B.y + A.y) + t * (B.y - A.y)) - mesh.ymin; // Shift si le maillage n'est pas centré

                        // 3. Remplissage de E
                        for (int i = 0; i < 3; ++i) {
                            for (int n = 0; n < N_modes; ++n) {

                                double cn_val = evaluate_c_1d(y_phys, h, n);

                                complex<double> val = phi_1d[i] * cn_val * w * detJac;
                                
                                E(nodes_global[i], n) += val;
                            }
                        }
                    }
                }
            }
        }
        return E;
    }

    // -------------------------------------------------------------------------
    // Assemblage matrice D
    // -------------------------------------------------------------------------

    void compute_D(FullMatrix<complexe>& D, int N_modes, double h, double k0){
            for (int j =0; j<N_modes; j++) {
                D(j,j) = complex<double>(0.0,1.0)* compute_beta(k0, h, j);
            }       
    }

    // -------------------------------------------------------------------------
    // Assemblage matrice T
    // -------------------------------------------------------------------------

    void T_matrix(ProfileMatrix<complexe>& K, FullMatrix<complexe>& E, FullMatrix<complexe>& D, double h, int boundary_tag, double factor = -1.0){

        int Ndof = E.rows();
        int Nmodes = D.rows();

        // Optimisation : Identifier les DOFs actifs (ceux sur le bord)
        // E est creuse (non nulle seulement sur le bord), on évite la boucle N^2
        vector<int> active_dofs;
        active_dofs.reserve(Ndof / 10); // Estimation
        for(int i=0; i<Ndof; ++i) {
            for(int n=0; n<Nmodes; ++n) {
                if(abs(E(i,n)) > 1e-14) { active_dofs.push_back(i); break; }
            }
        }

        for (int i : active_dofs) {
            
            for (int j : active_dofs) { 
                if (j > i) continue; // Symétrie : j <= i
                complexe val_T_ij = 0.0;
                
                for (int n = 0; n < Nmodes; ++n) {
                    // Calcul à la volée pour économiser la matrice temporaire ED
                    val_T_ij += E(i, n) * D(n, n) * E(j, n);
                }

                if (abs(val_T_ij) > 1e-14) {
                    K(i, j) += factor * val_T_ij;
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Assemblage vecteur G
    // -------------------------------------------------------------------------

    vector<complexe> assemble_source_vector(const MeshP2& mesh, 
                                               const FullMatrix<complexe>& E,
                                               int n_inc, 
                                               double k0, double L, double coef) {
    
        int Ndof = mesh.ndof();
        vector<complexe> G(Ndof, 0.0);
        
        complexe beta = compute_beta(k0, mesh.Ly, n_inc);
        complexe coeff = -2.0 * complex<double>(0,1) * beta * exp(coef*complex<double>(0,1) * beta * L);

        // Remplissage du vecteur

        for (int i = 0; i < Ndof; ++i) {
            G[i] = coeff * E(i, n_inc); 
        }
        return G;
    }
};

#endif