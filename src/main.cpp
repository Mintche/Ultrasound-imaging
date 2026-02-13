#include "fem.hpp"
#include "linear_sampling.hpp" // Ta classe développée précédemment
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

// DEFINE PARAMETERS
#define k0 30.0
#define kd (k0 * 3)
#define N_MODES 20
#define EPSILON 1e-6 // Paramètre de régularisation Tikhonov

// Tags physiques (à adapter selon le .msh)
#define TAG_LEFT 11
#define TAG_RIGHT 12
#define TAG_DEFECT 2

int main(int argc, char** argv) {

    if(argc < 2){
        printf("Usage: %s <mesh.msh>\n", argv[0]);
        return 1;
    }

    // -------------------------------------------------------------------------
    // INITIALISATION FEM & MATRICES SYSTEME
    // -------------------------------------------------------------------------
    printf("--- Initialisation FEM ---\n");
    MeshP2 mesh;
    mesh.read_msh_v2_ascii(argv[1], {TAG_DEFECT});

    double h = mesh.Ly;
    double L = mesh.Lx / 2.0;
    double x_source_gauche = mesh.xmin;
    double x_source_droite = mesh.xmax;
    int tag_left = TAG_LEFT;
    int tag_right = TAG_RIGHT;

    printf("H = %f | L = %f | Ndof : %zu\n", h, L, mesh.ndof());

    // Calcul du profil et matrices
    vector<size_t> profile = Fem::compute_profile_enhanced(mesh, {tag_left, tag_right});
    ProfileMatrix<complexe> K(profile);

    Fem::A_matrix(mesh, K, 1.0);                 // Matrice de Raideur
    Fem::B_matrix(mesh, K, k0, kd, -1.0);        // Matrice de Masse
    
    // Matrices de projection et diagonale D
    FullMatrix<complexe> E_minus = Fem::compute_E(mesh, N_MODES, tag_left, k0);
    FullMatrix<complexe> E_plus  = Fem::compute_E(mesh, N_MODES, tag_right, k0);
    FullMatrix<complexe> D(N_MODES, N_MODES);
    Fem::compute_D(D, N_MODES, h, k0);

    // Conditions DtN (Dirichlet-to-Neumann)
    Fem::T_matrix(K, E_minus, D, h, tag_left, -1.0);
    Fem::T_matrix(K, E_plus, D, h, tag_right, -1.0);

    // FACTORISATION DU SYSTEME FEM (On le fait une seule fois !)
    printf("Factorisation LDL^T de la matrice FEM...\n");
    K.factorize(); 

    // -------------------------------------------------------------------------
    // SIMULATION ET CONSTRUCTION DES MATRICES S (Scattering)
    // -------------------------------------------------------------------------
    printf("--- Simulation des ondes incidentes ---\n");
    
    // Matrices de Scattering S_Mesure_Source
    FullMatrix<complexe> S_LL(N_MODES, N_MODES); // Mesure Gauche, Source Gauche
    FullMatrix<complexe> S_RL(N_MODES, N_MODES); // Mesure Droite, Source Gauche
    FullMatrix<complexe> S_LR(N_MODES, N_MODES); // Mesure Gauche, Source Droite
    FullMatrix<complexe> S_RR(N_MODES, N_MODES); // Mesure Droite, Source Droite

    vector<complexe> U(mesh.ndof());
    vector<complexe> u_s(mesh.ndof());
    vector<complexe> proj_plus(N_MODES), proj_minus(N_MODES);

    // --- Cas 1 : Ondes incidentes depuis la Gauche (Source Plus / u_n^+) ---
    // Correspond aux colonnes de droite de F (F+-) et (F--)
    for (int n = 0; n < N_MODES; ++n) {
        vector<complexe> G = Fem::assemble_source_vector(mesh, E_minus, n, k0, x_source_gauche, 1.0);
        K.solve(U, G); 

        // Extraction champ diffracté (incident depuis gauche -> direction = 1.0)
        LinearSampling::compute_boundary_u_s(mesh, n, tag_left, tag_right, 1.0, x_source_gauche, k0, h, u_s, U);

        LinearSampling::compute_projection(u_s, E_plus, E_minus, proj_plus, proj_minus);

        for(int m=0; m<N_MODES; ++m) {
            S_RL(n, m) = proj_plus[m];  // (Un+)+ : Mesure Droite
            S_LL(n, m) = proj_minus[m]; // (Un+)- : Mesure Gauche
        }
    }

    // --- Cas 2 : Ondes incidentes depuis la Droite (Source Minus / u_n^-) ---
    // Correspond aux colonnes de gauche de F (F++) et (F-+)
    for (int n = 0; n < N_MODES; ++n) {
        vector<complexe> G = Fem::assemble_source_vector(mesh, E_plus, n, k0, x_source_droite, -1.0); 
        K.solve(U, G);

        // Extraction champ diffracté (incident depuis droite -> direction = -1.0)
        LinearSampling::compute_boundary_u_s(mesh, n, tag_left, tag_right, -1.0, x_source_droite, k0, h, u_s, U);

        LinearSampling::compute_projection(u_s, E_plus, E_minus, proj_plus, proj_minus);

        for(int m=0; m<N_MODES; ++m) {
            S_RR(n, m) = proj_plus[m];  // (Un-)+ : Mesure Droite
            S_LR(n, m) = proj_minus[m]; // (Un-)- : Mesure Gauche
        }
    }

    // -------------------------------------------------------------------------
    // CONSTRUCTION DU SYSTEME LINEAR SAMPLING
    // -------------------------------------------------------------------------
    printf("--- Construction LSM ---\n");

    FullMatrix<complexe> F = LinearSampling::compute_F(S_LL, S_RL, S_LR, S_RR, N_MODES, k0, h, L);

    // Construction de la matrice normale régularisée : M = F* F + epsilon * I

    FullMatrix<complexe> F_adj = F.adjoint();
    FullMatrix<complexe> I(2*N_MODES,2*N_MODES);
    for (int i = 0; i < 2*N_MODES;i++) I(i,i) = 1.0;
    FullMatrix<complexe> M = F_adj * F + complex(EPSILON,0.0)*I;

    M.factorize();

    // -------------------------------------------------------------------------
    // BOUCLE D'IMAGERIE
    // -------------------------------------------------------------------------
    printf("--- Calcul de l'image (Grille de points z) ---\n");

    // Paramètres de la grille d'imagerie
    int grid_nx = 100;
    int grid_ny = 50;
    // On scanne l'intérieur du maillage (marges de sécurité pour éviter les singularités aux bords)
    double x_scan_min = mesh.xmin + 0.05; double x_scan_max = mesh.xmax - 0.05;
    double y_scan_min = mesh.ymin + 0.05; double y_scan_max = mesh.ymax - 0.05;
    
    // Buffer pour stocker les résultats avant l'écriture fichier
    std::vector<double> indicators(grid_nx * grid_ny);

    printf("Calcul en cours sur %d points...\n", grid_nx * grid_ny);

        // Chaque thread a ses propres vecteurs de travail pour éviter les data-races
        vector<complexe> H_z(2*N_MODES);
        vector<complexe> F_adj_G(2*N_MODES);

        for(int i = 0; i < grid_nx; ++i) {
            for(int j = 0; j < grid_ny; ++j) {
                double z1 = x_scan_min + i * (x_scan_max - x_scan_min) / (grid_nx - 1);
                double z2 = y_scan_min + j * (y_scan_max - y_scan_min) / (grid_ny - 1);

                // 1. Calcul du second membre G_z
                vector<complexe> G_z = LinearSampling::assemble_Gz(mesh, N_MODES, z1, z2 - mesh.ymin, mesh.xmin, mesh.xmax, k0, h);

                // 2. Calcul du RHS : F* G_z
                F_adj_G = F_adj * G_z;

                // 3. Résolution (M est déjà factorisée)
                M.solve(H_z, F_adj_G);

                // 4. Calcul de la norme et de l'indicateur
                double norm_Hz = norm(H_z); 
                indicators[i * grid_ny + j] = 1.0 / (norm_Hz + 1e-15);
            }
        }

    // Écriture finale
    ofstream file("image_lsm.txt");
    file << grid_nx << " " << grid_ny << "\n";
    for(int i = 0; i < grid_nx; ++i) {
        double z1 = x_scan_min + i * (x_scan_max - x_scan_min) / (grid_nx - 1);
        for(int j = 0; j < grid_ny; ++j) {
            double z2 = y_scan_min + j * (y_scan_max - y_scan_min) / (grid_ny - 1);
            file << z1 << " " << z2 << " " << indicators[i * grid_ny + j] << "\n";
        }
    }
    file.close();

    printf("Termine. Resultats dans 'image_lsm.txt'.\n");
    return 0;
}