#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <string>
#include <cmath>

#include "fem.hpp"
#include "linear_sampling.hpp"

using namespace std;
using namespace usim;

int main(int argc, char** argv) {
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <mesh.msh> freqs" << endl;
        cout << "Exemple: " << argv[0] << " ../data/test_ultrasound_defaut_centre.msh f0 f1 f2..." << endl;
        return 1;
    }

    string mesh_file = argv[1];
    MeshP2 mesh;
    int tag_defect = 2;
    
    try {
        mesh.read_msh_v2_ascii(mesh_file, {tag_defect});
    } catch (const std::exception& e) {
        cerr << "Erreur de lecture du maillage: " << e.what() << endl;
        return 1;
    }
    
    Fem::reorder_mesh_rcm(mesh);

    int tag_left = 11;
    int tag_right = 12;
    int n_mode = 1;
    double c0 = 340;
    double contrast_ratio = 0.8;

    // Fichiers de sortie
    ofstream file_left("pinn_boundary_left_n0.csv");
    file_left << "f,k0,x,y,Re_U,Im_U\n";

    ofstream file_right("pinn_boundary_right_n0.csv");
    file_right << "f,k0,x,y,Re_U,Im_U\n";

    // Paramètres de fréquences (Curriculum Learning)
    int n_freqs = argc - 2;

    vector<size_t> profile = Fem::compute_profile_enhanced(mesh, {tag_left, tag_right});

    for (int i = 0; i < n_freqs; ++i) {
        
        double f = std::atof(argv[i+2]);
        if (f <= 0.0) {
            cerr << "Attention: Frequence invalide ignoree ('" << argv[i+2] << "')" << endl;
            continue;
        }
        
        double k0 = 2 * M_PI * f / c0;
        double kd = 2 * M_PI * f / (c0 * contrast_ratio);

        cout << "Generation PINN pour f = " << f << " Hz (k0 = " << k0 << ")..." << endl;

        ProfileMatrix<complexe> K(profile);
        Fem::A_matrix(mesh, K, 1.0);
        Fem::B_matrix(mesh, K, k0, kd, -1.0);

        // Même si on excite n=0, on a besoin de N_MODES suffisants pour absorber les modes évanescents
        int N_MODES = floor(mesh.Ly * k0 / M_PI) + 5; 

        FullMatrix<complexe> E_minus = Fem::compute_E(mesh, N_MODES, tag_left, k0);
        FullMatrix<complexe> E_plus  = Fem::compute_E(mesh, N_MODES, tag_right, k0);
        FullMatrix<complexe> D(N_MODES, N_MODES);
        Fem::compute_D(D, N_MODES, mesh.Ly, k0);

        Fem::T_matrix(K, E_minus, D, mesh.Ly, tag_left, -1.0);
        Fem::T_matrix(K, E_plus, D, mesh.Ly, tag_right, -1.0);

        K.factorize();

        // Excitation : Source Gauche, Mode n=0
        vector<complexe> G = Fem::assemble_source_vector(mesh, E_minus, n_mode, k0, mesh.xmin, 1.0);
        
        vector<complexe> U(mesh.ndof());
        K.solve(U, G);

        // Sauvegarde
        for (const auto& node : mesh.nodes) {
            if (node.ref == tag_left) {
                file_left << f << "," << k0 << "," << node.x << "," << node.y << "," << std::real(U[node.id]) << "," << std::imag(U[node.id]) << "\n";
            } else if (node.ref == tag_right) {
                file_right << f << "," << k0 << "," << node.x << "," << node.y << "," << std::real(U[node.id]) << "," << std::imag(U[node.id]) << "\n";
            }
        }
    }
    cout << "Donnees generees avec succes." << endl;
    return 0;
}