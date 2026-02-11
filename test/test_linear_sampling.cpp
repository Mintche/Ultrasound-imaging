#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "mesh.hpp"
#include "fem.hpp"
#include "linear_sampling.hpp"

using namespace std;
using namespace usim;

// Helper to check equality
bool check_close(complex<double> a, complex<double> b, double tol = 1e-9) {
    return abs(a - b) < tol;
}

int main() {
    cout << "--- Test Linear Sampling Module ---" << endl;

    // 1. Setup Mock Data
    int N_MODES = 5;
    double k0 = 10.0;
    double h = 1.0;
    double L = 2.0;

    // Create dummy S matrices (N x N)
    // S(incident, measure) convention based on main.cpp usage
    FullMatrix<complexe> S_pp(N_MODES, N_MODES);
    FullMatrix<complexe> S_pm(N_MODES, N_MODES);
    FullMatrix<complexe> S_mp(N_MODES, N_MODES);
    FullMatrix<complexe> S_mm(N_MODES, N_MODES);

    // Fill with known values to verify compute_F
    // We set S_mp to 1.0. In compute_F, F(m, n) uses S_m_p(n, m).
    // So if we set S_mp(n, m) = 1, then F should pick it up.
    S_mp.fill(complexe(1.0, 0.0));
    S_pp.fill(complexe(0.0, 0.0));
    S_mm.fill(complexe(0.0, 0.0));
    S_pm.fill(complexe(0.0, 0.0));

    // 2. Test compute_F
    cout << "Testing compute_F..." << endl;
    // Signature: compute_F(S_p_m, S_m_p, S_p_p, S_m_m, ...)
    FullMatrix<complexe> F = LinearSampling::compute_F(S_pm, S_mp, S_pp, S_mm, N_MODES, k0, h, L);

    // Check dimensions: F should be 2N x 2N
    if (F.rows() != (size_t)(2 * N_MODES)) {
         cout << "[FAIL] F dimensions incorrect. Expected " << 2*N_MODES << ", got " << F.rows() << endl;
         return 1;
    }

    // Check value for block F++ (top-left) which comes from S_mp
    // F(m, n) = C * S_m_p(n,m)
    // For n=0, m=0:
    complexe beta_0 = Fem::compute_beta(k0, h, 0);
    complexe exposant = complexe(0.,1.) * beta_0 * L;
    complexe constante = exp(exposant) / (complexe(0.,1.) * beta_0);
    
    complexe expected_F00 = constante * S_mp(0,0);
    
    if (check_close(F(0,0), expected_F00)) {
        cout << "[OK] F(0,0) calculation seems correct." << endl;
    } else {
        cout << "[FAIL] F(0,0) mismatch. Got " << F(0,0) << ", expected " << expected_F00 << endl;
    }

    // 3. Test assemble_Gz
    cout << "Testing assemble_Gz..." << endl;
    // assemble_Gz does not use mesh content in current implementation, so empty mesh is fine
    MeshP2 dummy_mesh;
    
    double z1 = 0.0; // center x
    double z2 = h/2.0; // center y
    
    vector<complexe> Gz = LinearSampling::assemble_Gz(dummy_mesh, N_MODES, z1, z2, L, k0, h);
    
    if (Gz.size() == (size_t)(2 * N_MODES)) {
        cout << "[OK] Gz size correct." << endl;
    } else {
        cout << "[FAIL] Gz size incorrect." << endl;
    }

    // 4. Test Full System Solve (Mock)
    cout << "Testing System Solve (M*g = RHS)..." << endl;
    // Construct M = F*F + eps*I
    FullMatrix<complexe> F_adj = F.adjoint();
    FullMatrix<complexe> I(2*N_MODES, 2*N_MODES);
    for(int i=0; i<2*N_MODES; ++i) I(i,i) = 1.0;
    
    double epsilon = 1e-4;
    FullMatrix<complexe> M = F_adj * F + complexe(epsilon, 0.0) * I;
    
    // Factorize
    M.factorize(); // LDLT
    cout << "[OK] Factorization successful." << endl;
    
    return 0;
}
