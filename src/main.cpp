#include "fem.hpp"
#include "linear_sampling.hpp"
#include <cstdio>
#include <iostream>
#include <vector>
#include <string>

// Paramètres physiques fixes (basés sur main_temp.cpp pour préserver la qualité de visualisation)
#define TAG_LEFT 11
#define TAG_RIGHT 12
#define TAG_DEFECT 2
#define BASE_K0 30.0
#define CONTRAST_RATIO 3.0 // kd = k0 * 3

int main(int argc, char** argv) {

    if(argc < 4){
        printf("Usage: %s <mesh.msh> <noise percentage> <n_frequencies>\n", argv[0]);
        printf("Exemple: %s mesh.msh 5 3\n", argv[0]);
        return 1;
    }

    std::string mesh_file = argv[1];
    double noise_percentage = atof(argv[2]);
    int n_freq = atoi(argv[3]);

    if (n_freq < 1) n_freq = 1;

    // -------------------------------------------------------------------------
    // INITIALISATION FEM
    // -------------------------------------------------------------------------
    printf("--- Initialisation FEM ---\n");
    MeshP2 mesh;
    try {
        mesh.read_msh_v2_ascii(mesh_file, {TAG_DEFECT});
    } catch (const std::exception& e) {
        std::cerr << "Erreur lors de la lecture du maillage: " << e.what() << std::endl;
        return 1;
    }
    
    // OPTIMISATION : Renumérotation RCM pour réduire la largeur de bande
    Fem::reorder_mesh_rcm(mesh);

    // Export pour vérification Matlab
    mesh.write_matlab_mesh_m("mesh_out.m");

    printf("Maillage charge : %zu noeuds, %zu triangles.\n", mesh.nodes.size(), mesh.triangles.size());
    printf("Parametres : Bruit = %.1f%%, Moyenne sur %d frequence(s).\n", noise_percentage, n_freq);

    // -------------------------------------------------------------------------
    // CALCUL LSM
    // -------------------------------------------------------------------------
    
    // Paramètres de la grille d'imagerie (identiques à main_temp)
    int grid_nx = 100;
    int grid_ny = 50;

    LinearSampling::compute_lsm_average(mesh, n_freq, BASE_K0, CONTRAST_RATIO, noise_percentage, 
                                        grid_nx, grid_ny, TAG_LEFT, TAG_RIGHT, "image_lsm.txt");

    return 0;
}