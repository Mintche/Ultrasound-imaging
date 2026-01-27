#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <cassert>

#include "fem.hpp"

using namespace std;
using namespace usim;

// Utilitaires de vérification
bool check_close(double val, double expected, double tol = 1e-9) {
    return std::abs(val - expected) < tol;
}

bool check_close(complexe val, complexe expected, double tol = 1e-9) {
    return std::abs(val - expected) < tol;
}

void test_shape_functions() {
    cout << "--- Test 1 : Fonctions de forme P2 ---" << endl;
    vector<double> phi(6);
    
    // Test au sommet 0 (0,0) -> phi[0] doit valoir 1
    Fem::evaluate_shape_functions(0.0, 0.0, phi);
    if(check_close(phi[0], 1.0) && check_close(phi[1], 0.0) && check_close(phi[3], 0.0))
        cout << "[OK] Valeurs au sommet 0 correctes." << endl;
    else
        cout << "[ECHEC] Valeurs au sommet 0 incorrectes." << endl;

    // Test au milieu de l'arête 0-1 (0.5, 0) -> phi[3] doit valoir 1
    Fem::evaluate_shape_functions(0.5, 0.0, phi);
    if(check_close(phi[3], 1.0) && check_close(phi[0], 0.0))
        cout << "[OK] Valeurs au noeud milieu correctes." << endl;
    else
        cout << "[ECHEC] Valeurs au noeud milieu incorrectes." << endl;

    // Partition de l'unité : somme(phi) = 1 partout
    Fem::evaluate_shape_functions(0.33, 0.33, phi);
    double sum = 0;
    for(double v : phi) sum += v;
    if(check_close(sum, 1.0))
        cout << "[OK] Partition de l'unite respectee." << endl;
    else
        cout << "[ECHEC] Partition de l'unite violee (somme=" << sum << ")." << endl;
}

void test_gradients() {
    cout << "\n--- Test 2 : Gradients P2 ---" << endl;
    FullMatrix<double> grads(6, 2);
    
    // La somme des gradients doit être nulle (grad(1) = 0)
    Fem::evaluate_gradients(0.25, 0.25, grads);
    double sum_dx = 0, sum_dy = 0;
    for(int i=0; i<6; ++i) {
        sum_dx += grads(i, 0);
        sum_dy += grads(i, 1);
    }
    
    if(check_close(sum_dx, 0.0) && check_close(sum_dy, 0.0))
        cout << "[OK] Somme des gradients nulle." << endl;
    else
        cout << "[ECHEC] Somme des gradients non nulle." << endl;
}

void test_assembly_single_triangle() {
    cout << "\n--- Test 3 : Assemblage (Rigidite & Masse) sur 1 Triangle ---" << endl;
    
    // Création manuelle d'un maillage avec 1 triangle (Triangle de référence)
    MeshP2 mesh;
    mesh.nodes.resize(6);
    // Sommets
    mesh.nodes[0] = {0.0, 0.0, 0};
    mesh.nodes[1] = {1.0, 0.0, 1};
    mesh.nodes[2] = {0.0, 1.0, 2};
    // Milieux
    mesh.nodes[3] = {0.5, 0.0, 3};
    mesh.nodes[4] = {0.5, 0.5, 4};
    mesh.nodes[5] = {0.0, 0.5, 5};
    
    TriangleP2 tri;
    tri.node_ids[0]=0; tri.node_ids[1]=1; tri.node_ids[2]=2;
    tri.node_ids[3]=3; tri.node_ids[4]=4; tri.node_ids[5]=5;
    tri.ref = 0;
    mesh.triangles.push_back(tri);
    
    // Calcul du profil
    vector<size_t> p = Fem::compute_profile(mesh);
    ProfileMatrix<complexe> K(p);
    ProfileMatrix<complexe> M(p);
    
    // Assemblage
    Fem::assemble_stiffness(mesh, K);
    Fem::assemble_mass(mesh, M, 1.0, 1.0); // k=1 partout
    
    // Vérification Rigidité : K * 1 = 0 (car grad(1) = 0)
    vector<complexe> ones(6, 1.0);
    vector<complexe> res_K = K * ones;
    bool k_ok = true;
    for(auto v : res_K) {
        if(abs(v) > 1e-9) k_ok = false;
    }
    if(k_ok) cout << "[OK] Matrice de Rigidite (noyau contient les constantes)." << endl;
    else cout << "[ECHEC] Matrice de Rigidite incorrecte." << endl;

    // Vérification Masse : 1^T * M * 1 = Int(1*1) = Aire = 0.5
    vector<complexe> res_M = M * ones;
    complexe total_mass = 0;
    for(auto v : res_M) total_mass += v;
    
    if(check_close(total_mass, 0.5)) 
        cout << "[OK] Matrice de Masse (integrale correcte)." << endl;
    else 
        cout << "[ECHEC] Matrice de Masse incorrecte (Masse=" << total_mass << ", attendu 0.5)." << endl;
}

void test_boundary_assembly() {
    cout << "\n--- Test 4 : Assemblage Bord (Matrice E) ---" << endl;
    
    // Maillage 1 triangle, l'arête 0 (noeuds 0-1) est sur le bord (tag 10)
    MeshP2 mesh;
    mesh.nodes.resize(6);
    mesh.nodes[0] = {0.0, 0.0, 0};
    mesh.nodes[1] = {1.0, 0.0, 1};
    mesh.nodes[2] = {0.0, 1.0, 2};
    mesh.nodes[3] = {0.5, 0.0, 3};
    mesh.nodes[4] = {0.5, 0.5, 4};
    mesh.nodes[5] = {0.0, 0.5, 5};
    
    TriangleP2 tri;
    tri.node_ids[0]=0; tri.node_ids[1]=1; tri.node_ids[2]=2;
    tri.node_ids[3]=3; tri.node_ids[4]=4; tri.node_ids[5]=5;
    
    // On marque l'arête 0 (entre noeuds 0 et 1) avec le tag 10
    tri.edge_ref[0] = 10; 
    tri.edge_ref[1] = 0; 
    tri.edge_ref[2] = 0;
    mesh.triangles.push_back(tri);
    
    vector<size_t> p = Fem::compute_profile(mesh);
    ProfileMatrix<complexe> E(p);
    
    // Assemblage sur le bord 10
    Fem::assemble_E(mesh, E, 10);
    
    // L'intégrale sur le bord (longueur 1) de la fonction constante 1 doit valoir 1.
    // Somme des éléments de E = 1^T * E * 1
    vector<complexe> ones(6, 1.0);
    vector<complexe> res = E * ones;
    complexe sum = 0;
    for(auto v : res) sum += v;
    
    if(check_close(sum, 1.0))
        cout << "[OK] Matrice de Bord (integrale longueur correcte)." << endl;
    else
        cout << "[ECHEC] Matrice de Bord incorrecte (Somme=" << sum << ", attendu 1.0)." << endl;
}

int main() {
    test_shape_functions();
    test_gradients();
    test_assembly_single_triangle();
    test_boundary_assembly();
    return 0;
}