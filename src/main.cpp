#include "fem.hpp"
#include "linear_sampling.hpp" // Ta classe développée précédemment
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

// DEFINE PARAMETERS
#define EPSILON 1e-6 // Paramètre de régularisation Tikhonov
#define NUM_K_POINTS 3
#define C 1500.0 // Vitesse du son dans le milieu de référence (m/s)
#define C_D 3000.0 // Vitesse du son dans le défaut (m/s)
// Tags physiques (à adapter selon le .msh)
#define TAG_LEFT 11
#define TAG_RIGHT 12
#define TAG_DEFECT 2


void compute_for_k(MeshP2 mesh,double k0, double kd, double h, double L, int tag_left, int tag_right, double x_source_gauche, double x_source_droite, double epsilon,
                                     int grid_nx, int grid_ny, vector<double>& indicators, double percentage, double noise_level) {
    int N_MODES = floor(h * min(k0,kd) / M_PI) + 5; 
    printf("Calcul pour k0 = %f (f = %f Hz) | N_modes = %d\n", k0, k0*C/(2*M_PI), N_MODES);                                    
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
    printf("Factorisation LDL^T de la matrice FEM \n");
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

        if (percentage != -1) LinearSampling::add_gaussian_noise(U, noise_level);

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

        if (percentage != -1) LinearSampling::add_gaussian_noise(U, noise_level);

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

   
    // On scanne l'intérieur du maillage (marges de sécurité pour éviter les singularités aux bords)
    double x_scan_min = mesh.xmin + 0.05; double x_scan_max = mesh.xmax - 0.05;
    double y_scan_min = mesh.ymin + 0.05; double y_scan_max = mesh.ymax - 0.05;
    


    printf("Calcul en cours sur %d points...\n", grid_nx * grid_ny);

    for(int i = 0; i < grid_nx; ++i) {
        for(int j = 0; j < grid_ny; ++j) {
                // Chaque thread a ses propres vecteurs de travail
                vector<complexe> H_z(2*N_MODES);
                vector<complexe> F_adj_G(2*N_MODES);

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

 }

void compute_average_image(MeshP2 mesh, int num_k_points, vector<double> &k_values, vector<double> &kd_values, double h, double L, int tag_left, int tag_right, double x_source_gauche, double x_source_droite, double epsilon,
                                     int grid_nx, int grid_ny,double percentage, double noise_level) {

    // Buffer pour stocker les images intermédiaires
    vector<vector<double>> all_indicators(num_k_points, vector<double>(grid_nx * grid_ny));

    // Boucle sur les différentes fréquences
    for (int idx = 0; idx < num_k_points; ++idx) {
        double k0 = k_values[idx];
        double kd = kd_values[idx];


        // Calcul des matrices et de l'image pour cette fréquence
        compute_for_k(mesh, k0, kd, h, L, tag_left, tag_right, x_source_gauche, x_source_droite, epsilon,
                     grid_nx, grid_ny, all_indicators[idx], percentage, noise_level);
    }

    // Moyenne des images intermédiaires
    vector<double> average_indicators(grid_nx * grid_ny, 0.0);
    for (int i = 0; i < grid_nx * grid_ny; ++i) {
        double sum = 0.0;
        for (int idx = 0; idx < num_k_points; ++idx) {
            sum += all_indicators[idx][i];
        }
        average_indicators[i] = sum / num_k_points;
    }

    // On scanne l'intérieur du maillage (marges de sécurité pour éviter les singularités aux bords)
    double x_scan_min = mesh.xmin + 0.05; double x_scan_max = mesh.xmax - 0.05;
    double y_scan_min = mesh.ymin + 0.05; double y_scan_max = mesh.ymax - 0.05;

    //Ecriture finale de l'image moyenne
    ofstream file("average_image_lsm.txt");
    file << grid_nx << " " << grid_ny << "\n";
    for(int i = 0; i < grid_nx; ++i) {
        double z1 = x_scan_min + i * (x_scan_max - x_scan_min) / (grid_nx - 1);
        for(int j = 0; j < grid_ny; ++j) {
            double z2 = y_scan_min + j * (y_scan_max - y_scan_min) / (grid_ny - 1);
            file << z1 << " " << z2 << " " << average_indicators[i * grid_ny + j] << " ";
        }
    }
    file.close();
    printf("Termine. Resultats dans 'average_image_lsm.txt'.\n");
}


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
    
    // OPTIMISATION : Renumérotation RCM pour réduire la largeur de bande
    Fem::reorder_mesh_rcm(mesh);

    mesh.write_matlab_mesh_m("mesh_out.m");
    //mesh.write_defect_coords_txt("defect_coords.txt"); si python

    int num_k_points = 3; // Nombre de points dans la plage de fréquences
    double h = mesh.Ly;
    double L = mesh.Lx / 2.0;
    double x_source_gauche = mesh.xmin;
    double x_source_droite = mesh.xmax;
    int tag_left = TAG_LEFT;
    int tag_right = TAG_RIGHT;
     // Paramètres de la grille d'imagerie
    int grid_nx = 100;
    int grid_ny = 50;
    double percentage = atof(argv[2]);
    double noise_level = (percentage == 0) ? -1 : percentage*(1/sqrt(2*h)); // Niveau de bruit


    
    // On génère une liste de k0 et kd pour une plage [k_max/10, K_max] avec k0=2*pi*f/c et kd=2*pi*f/c_d 
    //double k_max = 2.0 * M_PI / (h / 6.0); // k_max = 2*pi / (h/6)
    vector<double> kd_values(num_k_points);
    vector<double> k_values(num_k_points);
   /* double cmin=c;
    if (c_d < c) cmin = c_d;
    for (int i = 0; i < num_k_points; ++i) {
        double f = (i + 1) * k_max * cmin/ (2.0*M_PI * num_k_points); // Fréquence linéairement espacée
        k_values[i] = 2.0 * M_PI * f / c; // k0 = 2*pi*f/c
        kd_values[i] = 2.0 * M_PI * f / c_d; // kd = 2*pi*f/c_d
    }*/
    k_values[0] = 30;
    k_values[1] = 30;
    k_values[2] = 36;

    kd_values[0] = 3*30;
    kd_values[1] = 3*30;
    kd_values[2] = 3*36;


    printf("H = %f | L = %f | Ndof : %zu\n", h, L, mesh.ndof());

    compute_average_image(mesh, num_k_points,k_values,kd_values, h, L, tag_left, tag_right, x_source_gauche , x_source_droite, EPSILON,
                     grid_nx, grid_ny, percentage, noise_level);
    return 0;
}
