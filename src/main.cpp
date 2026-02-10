#include "fem.hpp"
#include "linear_sampling.hpp" // Ta classe développée précédemment
#include <cstdio>
#include <fstream>
#include <iostream>

// DEFINE PARAMETERS
#define k0 30.0
#define kd (k0 * 3)
#define N_MODES 20
#define EPSILON 1e-4 // Paramètre de régularisation Tikhonov

int main(int argc, char** argv) {

    if(argc < 2){
        printf("Usage: %s <mesh.msh>\n", argv[0]);
        return 1;
    }

    // -------------------------------------------------------------------------
    // 1. INITIALISATION FEM & MATRICES SYSTEME
    // -------------------------------------------------------------------------
    printf("--- Initialisation FEM ---\n");
    MeshP2 mesh;
    mesh.read_msh_v2_ascii(argv[1], {2});

    double h = mesh.Ly;
    double L = mesh.Lx / 2.0;
    int tag_left = 11;
    int tag_right = 12;

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
    Fem::compute_D(D, 0, N_MODES, h, k0);

    // Conditions DtN (Dirichlet-to-Neumann)
    Fem::T_matrix(K, E_minus, D, h, tag_left, -1.0);
    Fem::T_matrix(K, E_plus, D, h, tag_right, -1.0);

    // FACTORISATION DU SYSTEME FEM (On le fait une seule fois !)
    printf("Factorisation LDL^T de la matrice FEM...\n");
    K.factorize(); 

    // -------------------------------------------------------------------------
    // 2. SIMULATION ET CONSTRUCTION DES MATRICES S (Scattering)
    // -------------------------------------------------------------------------
    printf("--- Simulation des ondes incidentes ---\n");
    
    // Matrices de Scattering S(mesure, incident)
    // S_m_p : Mesure Minus (Gauche), Incident Plus (depuis Gauche vers Droite)
    FullMatrix<complexe> S_m_p(N_MODES, N_MODES); 
    FullMatrix<complexe> S_p_p(N_MODES, N_MODES);
    FullMatrix<complexe> S_m_m(N_MODES, N_MODES);
    FullMatrix<complexe> S_p_m(N_MODES, N_MODES);

    vector<complexe> U(mesh.ndof());
    vector<complexe> u_s(mesh.ndof());
    vector<complexe> proj_plus(N_MODES), proj_minus(N_MODES);

    // --- Cas 1 : Ondes incidentes depuis la Gauche (Phi+) ---
    for (int n = 0; n < N_MODES; ++n) {
        // Source pour le mode n entrant par la gauche
        vector<complexe> G = Fem::assemble_source_vector(mesh, E_minus, n, k0, L, -1.0);
        
        // Résolution (utilise la factorisation pré-calculée)
        K.solve(U, G); 

        // Extraction champ diffracté u^s = U_tot - Phi_inc
        // true pour incident_from_left
        LinearSampling::compute_boundary_u_s(mesh, n, tag_left, tag_right, 1.0, L, k0, h, u_s, U);

        // Projections sur les modes
        LinearSampling::compute_projection(u_s, E_plus, E_minus, proj_plus, proj_minus);

        // Remplissage des colonnes n des matrices S
        for(int m=0; m<N_MODES; ++m) {
            S_p_p(m, n) = proj_plus[m];  // Mesure Droite, Inc Gauche
            S_p_m(m, n) = proj_minus[m]; // Mesure Gauche, Inc Gauche
        }
    }

    // --- Cas 2 : Ondes incidentes depuis la Droite (Phi-) ---
    for (int n = 0; n < N_MODES; ++n) {
        // Source pour le mode n entrant par la droite (utilise E_plus)
        vector<complexe> G = Fem::assemble_source_vector(mesh, E_plus, n, k0, L, 1.0); 
        
        K.solve(U, G);

        // Extraction champ diffracté (false pour incident_from_right)
        LinearSampling::compute_boundary_u_s(mesh, n, tag_left, tag_right, -1.0, L, k0, h, u_s, U);

        LinearSampling::compute_projection(u_s, E_plus, E_minus, proj_plus, proj_minus);

        for(int m=0; m<N_MODES; ++m) {
            S_p_m(m, n) = proj_plus[m];
            S_m_m(m, n) = proj_minus[m];
        }
    }

    // -------------------------------------------------------------------------
    // 3. CONSTRUCTION DU SYSTEME LINEAR SAMPLING
    // -------------------------------------------------------------------------
    printf("--- Construction LSM ---\n");

    // Construction de F (2N x 2N)
    FullMatrix<complexe> F = LinearSampling::compute_F(S_m_p, S_p_p, S_m_m, S_p_m, N_MODES, k0, h, L);

    // Construction de la matrice normale régularisée : M = F* F + epsilon * I
    // On suppose que ta classe FullMatrix gère l'adjoint et le produit
    FullMatrix<complexe> F_adj = F.adjoint();
    FullMatrix<complexe> I(2*N_MODES,2*N_MODES);
    for (int i = 0; i < 2*N_MODES;i++) I(i,i) = 1.0;
    FullMatrix<complexe> M = F_adj * F + complex(EPSILON,0.0)*I;

    M.factorize();

    // -------------------------------------------------------------------------
    // 4. BOUCLE D'IMAGERIE
    // -------------------------------------------------------------------------
    printf("--- Calcul de l'image (Grille de points z) ---\n");

    // Paramètres de la grille d'imagerie
    int grid_nx = 100;
    int grid_ny = 50;
    double x_min = -L + 0.1, x_max = L - 0.1;
    double y_min = 0.0 + 0.05, y_max = h - 0.05;

    ofstream file("image_lsm.txt");
    file << grid_nx << " " << grid_ny << "\n"; // Header

    vector<complexe> H_z(2*N_MODES);
    vector<complexe> RHS(2*N_MODES);

    for(int i = 0; i < grid_nx; ++i) {
        double z1 = x_min + i * (x_max - x_min) / (grid_nx - 1);
        
        for(int j = 0; j < grid_ny; ++j) {
            double z2 = y_min + j * (y_max - y_min) / (grid_ny - 1);

            // 1. Calcul du second membre G_z (Vecteur Green théorique)
            vector<complexe> G_z = LinearSampling::assemble_Gz(mesh, N_MODES, z1, z2, L, k0, h);

            // 2. Calcul du RHS du système normal : F* G_z
            
            vector<complexe> F_adj_G = F_adj * G_z;

            // 3. Résolution du système (M * H_z = RHS)
            M.solve(H_z, F_adj_G);

            // 4. Calcul de la norme
            double norm_Hz = 0.0;
            for(const auto& val : H_z) {
                norm_Hz += norm(val); // norm renvoie |z|^2
            }
            norm_Hz = sqrt(norm_Hz);

            // 5. Stockage (Indicateur = 1 / ||Hz||)
            double indicator = 1.0 / (norm_Hz + 1e-15); // Eviter div par 0
            
            file << z1 << " " << z2 << " " << indicator << "\n";
        }
    }
    file.close();

    printf("Termine. Resultats dans 'image_lsm.txt'.\n");
    return 0;
}