#include "fem.hpp"
#include <cstdio>

//DEF

#define k0 1.0
#define kd (k0*1.5)
#define N 20

int main(int argc,char** argv){

    //import du maillage, creation des matrice, et de la solution.

    printf("Import du maillage...\n");

    MeshP2 mesh;

    mesh.read_msh_v2_ascii("data/test_ultrasound.msh",{2});

    double h = mesh.Ly;

    printf("Hauteur du guide (h) : %f | Ndof : %zu\n", h, mesh.ndof());

    vector<size_t> profile = Fem::compute_profile_enhanced(mesh,{11, 12}); // 11 bord gauche 12 bord droite

    ProfileMatrix<complexe> K(profile);

    Fem::A_matrix(mesh,K,1.0); //on ajoute A

    Fem::B_matrix(mesh,K,k0,kd,-1.0); //on soustrait B

    // --- Conditions aux limites transparentes (DtN) ---
    // On calcule les termes de bord pour simuler un guide infini
    int N_modes = N;
    int tag_left = 11;  // Tag du bord gauche dans le .msh
    int tag_right = 12; // Tag du bord droit dans le .msh

    // Matrices de projection E (liens entre noeuds du bord et modes)
    FullMatrix<complexe> E_minus = Fem::compute_E(mesh, N_modes, tag_left, k0);
    FullMatrix<complexe> E_plus  = Fem::compute_E(mesh, N_modes, tag_right, k0);

    // Matrice diagonale D (contient les coefficients de propagation i*beta_n)
    FullMatrix<complexe> D(N_modes, N_modes);
    Fem::compute_D(D, 0, N_modes, h, k0);

    // Ajout des matrices T = E * D * E^T à la matrice globale K
    // Le facteur est -1.0 car le terme de bord passe à droite ou est soustrait dans la formulation faible
    Fem::T_matrix(K, E_minus, D, h, tag_left, -1.0);
    Fem::T_matrix(K, E_plus, D, h, tag_right, -1.0);

    // --- Assemblage du second membre (Source) ---
    // On envoie le mode fondamental (n=0) depuis la gauche
    double L = mesh.Lx / 2.0; // Demi-longueur du domaine
    vector<complexe> G = Fem::assemble_source_vector(mesh, E_minus, 0, k0, L);

    // --- Résolution du système KU = G ---
    printf("Resolution du systeme (Ndof=%zu)...\n", mesh.ndof());
    vector<complexe> U(mesh.ndof());
    K.solve(U, G); // Factorisation LDL^T + Descente/Remontée

    // --- Export pour Matlab ---
    printf("Export des resultats...\n");
    mesh.write_matlab_mesh_m("mesh_out.m");       // Exporte la géométrie
    mesh.write_matlab_field_m("solution_out.m", U); // Exporte le champ U

    printf("Termine.\n");
    return 0;
}