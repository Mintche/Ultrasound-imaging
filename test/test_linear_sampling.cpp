#include "mesh.hpp"
#include "linear_sampling.hpp"

// -----------------------------------------------------------------------------
// Minimal smoke test for the imaging module.
//
// This program does NOT depend on the FEM part. It only demonstrates the
// expected data flow once you have scattered fields from your simulator.
//
// Usage (from repo root):
//   g++ -O2 -std=c++17 -Iinclude test/test_linear_sampling.cpp -o test_lsm
//   ./test_lsm
// -----------------------------------------------------------------------------

using namespace usim;

int main(){
    // 1) Read mesh and infer geometric parameters
    MeshP2 mesh;
    mesh.read_msh_v2_ascii("test/test_ultrasound.msh", {2});
    mesh.compute_bbox_and_dims();

    // In the statement, h is the waveguide height.
    const double h = mesh.Ly;
    // In the statement, the domain is centered x in [-L, +L].
    // Here the test mesh is in [0, Lx], so *centered coordinate* is x - Lx/2.
    const double L = 0.5 * mesh.Lx;

    // 2) Choose physics (example numbers)
    const double k0 = 20.0; // w/c0
    const int N = 5;        // modes 0..N

    // Boundary physical tags in test_ultrasound.geo:
    //  - 11 (left), 12 (right)
    const int tag_left  = 11;
    const int tag_right = 12;

    // 3) Build E matrices on the two measurement lines
    // Here we interpret:
    //   Sigma_- = left boundary, Sigma_+ = right boundary.
    const auto E_plus  = build_E_on_boundary(mesh, tag_right, N, h);
    const auto E_minus = build_E_on_boundary(mesh, tag_left,  N, h);

    // 4) Dummy scattered fields (replace by your FEM results)
    const int nm = N + 1;
    std::vector<std::vector<complexe>> Uplus_fields(nm), Uminus_fields(nm);
    for(int n=0;n<nm;++n){
        Uplus_fields[n]  = std::vector<complexe>(mesh.ndof(), complexe(0.0));
        Uminus_fields[n] = std::vector<complexe>(mesh.ndof(), complexe(0.0));
    }

    // 5) Convert nodal fields -> modal measurements
    const auto meas = build_measurements_from_fields(N, E_plus, E_minus, Uplus_fields, Uminus_fields);

    // 6) Build F and run an imaging grid
    LinearSampling lsm(N, k0, h, L, LinearSampling::Mode::FullScattering);
    lsm.build_F(meas);

    // We sample in centered coords: x in [-L, L], y in [0, h]
    lsm.image_grid(-L, L, 0.0, h, 60, 40, 1e-3, "lsm_image.xyz");

    return 0;
}
