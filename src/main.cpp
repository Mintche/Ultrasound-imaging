#include "fem.hpp"
#include <cstdio>

//DEFINE

#define k0 30.0
#define kd (k0*1.5)
#define N 20

int main(int argc,char** argv){

    //import du maillage, creation des matrice, et de la solution.

    printf("Import du maillage...\n");

    MeshP2 mesh;

    mesh.read_msh_v2_ascii(argv[1],{2});

    double h = mesh.Ly;
    int N_modes = N;
    int tag_left = 11;
    int tag_right = 12;

    printf("Hauteur du guide (h) : %f | Ndof : %zu\n", h, mesh.ndof());

    vector<size_t> profile = Fem::compute_profile_enhanced(mesh,{tag_left, tag_right});

    ProfileMatrix<complexe> K(profile);

    Fem::A_matrix(mesh,K,1.0); //on ajoute A

    Fem::B_matrix(mesh,K,k0,kd,-1.0); //on soustrait B

    // Matrices de projection E

    FullMatrix<complexe> E_minus = Fem::compute_E(mesh, N_modes, tag_left, k0);
    FullMatrix<complexe> E_plus  = Fem::compute_E(mesh, N_modes, tag_right, k0);

    // Matrice diagonale D (contient les coefficients de propagation i*beta_n)

    FullMatrix<complexe> D(N_modes, N_modes);
    Fem::compute_D(D, 0, N_modes, h, k0);

    // Ajout des matrices T

    Fem::T_matrix(K, E_minus, D, h, tag_left, -1.0);
    Fem::T_matrix(K, E_plus, D, h, tag_right, -1.0);

    // --- Assemblage du second membre (Source) ---
    // On envoie le mode n depuis la gauche vers la droite

    int n = 0;

    double L = mesh.Lx / 2.0; // Demi-longueur du domaine

    vector<complexe> G = Fem::assemble_source_vector(mesh, E_minus, n, k0, L);

    // --- Résolution du système KU = G ---
    printf("Resolution du systeme (Ndof=%zu)...\n", mesh.ndof());

    vector<complexe> U(mesh.ndof());

    K.solve(U, G);

    // --- Export pour Matlab ---
    printf("Export des resultats...\n");
    mesh.write_matlab_mesh_m("mesh_out.m");       // Exporte la géométrie
    mesh.write_matlab_field_m("solution_out.m", U); // Exporte le champ U

    printf("Termine.\n");
    return 0;
}