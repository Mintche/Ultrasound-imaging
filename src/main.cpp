#include "fem.hpp"
#include "linear_sampling.hpp"
#include <cstdio>
#include <iostream>
#include <vector>
#include <string>

constexpr int TAG_LEFT = 11;
constexpr int TAG_RIGHT = 12;
constexpr int TAG_DEFECT = 2;
constexpr double BASE_K0 = 30.0;
constexpr double CONTRAST_RATIO = 3.0;

int main(int argc, char** argv) {

    if(argc < 4){
        std::cout << "Usage: "<< argv[0]<< "<mesh.msh> <noise percentage> <n_frequencies> "<< std::endl;
        std::cout <<"Exemple: "<< argv[0]<< " mesh.msh 5 3"<< std::endl;
        return 1;
    }

    std::string mesh_file = argv[1];
    double noise_percentage = atof(argv[2]);
    int n_freq = atoi(argv[3]);

    if (n_freq < 1) n_freq = 1;

    // -------------------------------------------------------------------------
    // INITIALISATION FEM
    // -------------------------------------------------------------------------
    std::cout << "--- Initialisation FEM ---\n";
    usim::MeshP2 mesh;
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

    std::cout << "Maillage charge : " <<  mesh.nodes.size() <<  " noeuds, " << mesh.nodes.size() << " triangles." <<std::endl;
    std::cout << "Parametres : Bruit = " << 100*noise_percentage<< " Moyenne sur " << n_freq << "frequence(s)";

    // -------------------------------------------------------------------------
    // CALCUL LSM
    // -------------------------------------------------------------------------

    LinearSampling::PhysicalParameters phys_params = {340, 15*1e2, 2*1e1, n_freq}; //c0 freq start frestep n_freq
    
    // Paramètres de la grille d'imagerie
    int grid_nx = 100;
    int grid_ny = 50;

    //LinearSampling::compute_lsm_average(mesh, n_freq, BASE_K0, CONTRAST_RATIO, noise_percentage, grid_nx, grid_ny, TAG_LEFT, TAG_RIGHT, "image_lsm.txt");

    LinearSampling::compute_lsm_physical(mesh,phys_params,5,noise_percentage, grid_nx, grid_ny, TAG_LEFT, TAG_RIGHT, "image_lsm.txt");

    return 0;
}