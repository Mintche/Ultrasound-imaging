#ifndef LINEAR_SAMPLING_HPP
#define LINEAR_SAMPLING_HPP

#include <cmath>
#include <complex>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "math.hpp" // FullMatrix, complexe
#include "mesh.hpp" // MeshP2

// -----------------------------------------------------------------------------
// Linear Sampling Method (LSM) - Imaging part
// -----------------------------------------------------------------------------
//
// This header implements the imaging algorithm described in the project notes:
//   - build the near-field operator matrix F from modal measurements,
//   - solve the Tikhonov-regularized system for many sampling points z,
//   - output an "image" value: log(1 / ||H_eps(z)||).
//
// It is designed to plug into your simulator once you can generate scattered
// fields for incident waveguide modes.
//
// Conventions used here follow section 2.2 of the statement:
//   - modes are indexed n=0..N
//   - beta_n = sqrt(k0^2 - (n*pi/h)^2) with branch Re>=0, Im>=0
//   - c_n(x2) are L2(0,h)-normalized cosines
//   - F is assembled from modal coefficients (U^±_n)^±_m of scattered fields
//
// IMPORTANT (geometry): The statement assumes a centered domain x1 in [-L,+L]
// with measurement lines Sigma_{+} = {+L}x(0,h) and Sigma_{-} = {-L}x(0,h).
// If your mesh is not centered (e.g. x in [0,Lx]), you can either:
//   (a) shift coordinates before feeding z1 into build_rhs_G(...), or
//   (b) choose L as half-width and express z1 in centered coordinates.
// -----------------------------------------------------------------------------

namespace usim {

// --------- small helpers ---------
inline complexe I_unit(){ return complexe(0.0, 1.0); }

inline double pi(){ return 3.141592653589793238462643383279502884; }

// L2-normalized cosine mode c_n(x2) on (0,h)
inline double c_n(int n, double x2, double h){
    if(h <= 0.0) throw std::invalid_argument("c_n: h must be > 0");
    const double a = (n==0) ? (1.0/std::sqrt(h)) : (std::sqrt(2.0/h));
    return a * std::cos(n * pi() * x2 / h);
}

// beta_n with branch Re>=0, Im>=0.
inline complexe beta_n(int n, double k0, double h){
    const double alpha = n * pi() / h;
    const complexe val = complexe(k0*k0 - alpha*alpha, 0.0);
    complexe b = std::sqrt(val);
    // ensure Re>=0
    if(std::real(b) < 0) b = -b;
    // ensure Im>=0
    if(std::imag(b) < 0) b = complexe(std::real(b), -std::imag(b));
    return b;
}

// --------- boundary integration for E matrix ---------

namespace detail_lsm {

struct EdgeKey {
    int a=0,b=0;
    EdgeKey() = default;
    EdgeKey(int i,int j){
        if(i<j){a=i;b=j;} else {a=j;b=i;}
    }
    bool operator==(const EdgeKey& o) const { return a==o.a && b==o.b; }
};

struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& e) const noexcept {
        return (static_cast<std::size_t>(e.a) << 32) ^ static_cast<std::size_t>(e.b);
    }
};

// 1D P2 Lagrange basis on [-1,1] with nodes at s=-1,0,1.
inline void p2_basis_1d(double s, double& phi_m1, double& phi_0, double& phi_p1){
    // l1(s) at -1, l2(s) at 0, l3(s) at +1
    phi_m1 = 0.5 * s * (s - 1.0);
    phi_0  = 1.0 - s*s;
    phi_p1 = 0.5 * s * (s + 1.0);
}

} // namespace detail_lsm

// Build E matrix (Ndof x (N+1)) for a given boundary physical tag.
//
// E(I,n) = \int_{Sigma} w_I(x) c_n(x2) ds(x)
//
// Here w_I restricted to an edge is a 1D quadratic Lagrange basis.
// We integrate each boundary edge with a 3-pt Gauss-Legendre rule on [-1,1].
inline FullMatrix<complexe> build_E_on_boundary(const MeshP2& mesh,
                                                int boundary_tag,
                                                int N,
                                                double h){
    const int ndof = static_cast<int>(mesh.ndof());
    const int nm   = N + 1;
    FullMatrix<complexe> E(ndof, nm);

    // Gauss-Legendre 3-pt on [-1,1]
    const double s_q[3] = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
    const double w_q[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

    // Avoid double-counting if a boundary edge appears twice (shouldn't, but safe)
    std::unordered_set<detail_lsm::EdgeKey, detail_lsm::EdgeKeyHash> seen;

    for(const auto& T : mesh.triangles){
        // edges: 0->(1,2) mid=4, 1->(2,3) mid=5, 2->(3,1) mid=6
        for(int e=0;e<3;++e){
            if(T.edge_ref[e] != boundary_tag) continue;

            int v1=-1, v2=-1, vm=-1;
            if(e==0){ v1=T.node_ids[0]; v2=T.node_ids[1]; vm=T.node_ids[3]; }
            if(e==1){ v1=T.node_ids[1]; v2=T.node_ids[2]; vm=T.node_ids[4]; }
            if(e==2){ v1=T.node_ids[2]; v2=T.node_ids[0]; vm=T.node_ids[5]; }

            detail_lsm::EdgeKey key(v1,v2);
            if(seen.find(key)!=seen.end()) continue;
            seen.insert(key);

            const auto& P1 = mesh.nodes[static_cast<size_t>(v1)];
            const auto& P2 = mesh.nodes[static_cast<size_t>(v2)];
            const double dx = P2.x - P1.x;
            const double dy = P2.y - P1.y;
            const double len = std::sqrt(dx*dx + dy*dy);
            const double jac = 0.5 * len; // |dX/ds| on [-1,1]

            for(int q=0;q<3;++q){
                const double s = s_q[q];
                const double w = w_q[q];
                const double x2 = 0.5*(1.0 - s)*P1.y + 0.5*(1.0 + s)*P2.y;

                double phi_m1, phi_0, phi_p1;
                detail_lsm::p2_basis_1d(s, phi_m1, phi_0, phi_p1);

                for(int n=0;n<nm;++n){
                    const double cn = c_n(n, x2, h);
                    const complexe f = complexe(cn * w * jac, 0.0);
                    // node at s=-1 -> v1 ; s=0 -> vm ; s=+1 -> v2
                    E(v1,n) += f * phi_m1;
                    E(vm,n) += f * phi_0;
                    E(v2,n) += f * phi_p1;
                }
            }
        }
    }

    return E;
}

// Modal coefficients on one boundary: coeff[n] = (u|Sigma, c_n) ≈ Σ_I U_I * E(I,n)
inline std::vector<complexe> modal_decomposition(const FullMatrix<complexe>& E,
                                                 const std::vector<complexe>& U,
                                                 int nm){
    const int ndof = static_cast<int>(U.size());
    std::vector<complexe> a(static_cast<size_t>(nm), complexe(0.0));
    for(int n=0;n<nm;++n){
        complexe s(0.0);
        for(int I=0; I<ndof; ++I) s += U[static_cast<size_t>(I)] * E(I,n);
        a[static_cast<size_t>(n)] = s;
    }
    return a;
}

// Container for the modal measurements needed to build F.
// Each matrix is (nm x nm): rows m, cols n.
struct ModalMeasurements {
    int N = 0; // highest mode index
    // (U^-_n)^+_m, (U^+_n)^+_m, (U^-_n)^-_m, (U^+_n)^-_m
    FullMatrix<complexe> Uminus_on_plus;
    FullMatrix<complexe> Uplus_on_plus;
    FullMatrix<complexe> Uminus_on_minus;
    FullMatrix<complexe> Uplus_on_minus;

    explicit ModalMeasurements(int N_)
        : N(N_),
          Uminus_on_plus(N_+1, N_+1),
          Uplus_on_plus(N_+1, N_+1),
          Uminus_on_minus(N_+1, N_+1),
          Uplus_on_minus(N_+1, N_+1) {}
};

// Build modal measurements from nodal scattered fields.
//
// Inputs:
//   - E_plus : E matrix on Sigma_{+} (x=+L)  [ndof x nm]
//   - E_minus: E matrix on Sigma_{-} (x=-L)  [ndof x nm]
//   - Uplus_fields[n]  : scattered field for incident Phi^+_n (size ndof)
//   - Uminus_fields[n] : scattered field for incident Phi^-_n (size ndof)
//
// Outputs matrices have rows m (observed mode), cols n (incident mode).
inline ModalMeasurements build_measurements_from_fields(
        int N,
        const FullMatrix<complexe>& E_plus,
        const FullMatrix<complexe>& E_minus,
        const std::vector<std::vector<complexe>>& Uplus_fields,
        const std::vector<std::vector<complexe>>& Uminus_fields){

    const int nm = N + 1;
    if(static_cast<int>(Uplus_fields.size()) != nm || static_cast<int>(Uminus_fields.size()) != nm)
        throw std::invalid_argument("build_measurements_from_fields: need (N+1) fields for + and - incidents");

    ModalMeasurements meas(N);

    for(int n=0;n<nm;++n){
        const auto a_plus_on_plus  = modal_decomposition(E_plus,  Uplus_fields[static_cast<size_t>(n)],  nm);
        const auto a_plus_on_minus = modal_decomposition(E_minus, Uplus_fields[static_cast<size_t>(n)],  nm);
        const auto a_minus_on_plus  = modal_decomposition(E_plus,  Uminus_fields[static_cast<size_t>(n)], nm);
        const auto a_minus_on_minus = modal_decomposition(E_minus, Uminus_fields[static_cast<size_t>(n)], nm);

        for(int m=0;m<nm;++m){
            meas.Uplus_on_plus(m,n)   = a_plus_on_plus[static_cast<size_t>(m)];
            meas.Uplus_on_minus(m,n)  = a_plus_on_minus[static_cast<size_t>(m)];
            meas.Uminus_on_plus(m,n)  = a_minus_on_plus[static_cast<size_t>(m)];
            meas.Uminus_on_minus(m,n) = a_minus_on_minus[static_cast<size_t>(m)];
        }
    }

    return meas;
}

// --------- LSM core ---------

class LinearSampling {
public:
    enum class Mode { FullScattering, BackScattering };

    LinearSampling(int N_, double k0_, double h_, double L_, Mode mode_ = Mode::FullScattering)
        : N(N_), nm(N_+1), k0(k0_), h(h_), L(L_), mode(mode_),
          betas(static_cast<size_t>(nm), complexe(0.0)),
          F( (mode_==Mode::FullScattering)? 2*nm : nm,
             (mode_==Mode::FullScattering)? 2*nm : nm )
    {
        if(N < 0) throw std::invalid_argument("LinearSampling: N must be >= 0");
        if(k0 <= 0.0) throw std::invalid_argument("LinearSampling: k0 must be > 0");
        if(h <= 0.0) throw std::invalid_argument("LinearSampling: h must be > 0");
        if(L <= 0.0) throw std::invalid_argument("LinearSampling: L must be > 0");

        for(int n=0;n<nm;++n) betas[static_cast<size_t>(n)] = beta_n(n, k0, h);
    }

    // Build F from measurements (eq. 2.2 of statement).
    // For BackScattering, only the F++ block is used.
    void build_F(const ModalMeasurements& meas){
        if(meas.N != N) throw std::invalid_argument("build_F: meas.N mismatch");

        const complexe ii = I_unit();
        if(mode == Mode::BackScattering){
            // F := F++ (nm x nm)
            for(int m=0;m<nm;++m){
                for(int n=0;n<nm;++n){
                    const complexe b = betas[static_cast<size_t>(n)];
                    F(m,n) = std::exp(ii*b*L) / (ii*b) * meas.Uminus_on_plus(m,n);
                }
            }
            return;
        }

        // Full scattering: F is 2nm x 2nm in blocks
        // [F++ F+-; F-+ F--]
        for(int m=0;m<nm;++m){
            for(int n=0;n<nm;++n){
                const complexe b = betas[static_cast<size_t>(n)];
                const complexe s = std::exp(ii*b*L) / (ii*b);

                const complexe Fpp = s * meas.Uminus_on_plus(m,n);
                const complexe Fpm = s * meas.Uplus_on_plus(m,n);
                const complexe Fmp = s * meas.Uminus_on_minus(m,n);
                const complexe Fmm = s * meas.Uplus_on_minus(m,n);

                // row/col indexing
                const int rp = m;
                const int rm = nm + m;
                const int cp = n;
                const int cm = nm + n;

                F(rp,cp) = Fpp;
                F(rp,cm) = Fpm;
                F(rm,cp) = Fmp;
                F(rm,cm) = Fmm;
            }
        }
    }

    // Compute the LSM indicator on a regular grid and write a 3-column file: x y val.
    // val = log(1 / ||H_eps(z)||).
    //
    // x and y are in the *centered* coordinate system of the statement.
    void image_grid(double x_min, double x_max,
                    double y_min, double y_max,
                    int nx, int ny,
                    double eps_reg,
                    const std::string& out_xyz_file) const
    {
        if(nx <= 1 || ny <= 1) throw std::invalid_argument("image_grid: nx,ny must be > 1");
        if(eps_reg <= 0.0) throw std::invalid_argument("image_grid: eps_reg must be > 0");

        const int dim = (mode==Mode::FullScattering) ? 2*nm : nm;

        // Build A = F^H F + eps I once
        FullMatrix<complexe> A(dim, dim);
        for(int i=0;i<dim;++i){
            for(int j=0;j<dim;++j){
                complexe s(0.0);
                for(int k=0;k<dim;++k){
                    s += std::conj(F(k,i)) * F(k,j);
                }
                if(i==j) s += eps_reg;
                A(i,j) = s;
            }
        }

        std::ofstream out(out_xyz_file);
        if(!out) throw std::runtime_error("image_grid: cannot write to " + out_xyz_file);
        out << "# x y val\n";

        const double dx = (x_max - x_min) / (nx - 1);
        const double dy = (y_max - y_min) / (ny - 1);

        for(int iy=0; iy<ny; ++iy){
            const double z2 = y_min + iy*dy;
            for(int ix=0; ix<nx; ++ix){
                const double z1 = x_min + ix*dx;
                const auto Gz = build_rhs_G(z1, z2);

                // b = F^H G
                std::vector<complexe> b(dim, complexe(0.0));
                for(int i=0;i<dim;++i){
                    complexe s(0.0);
                    for(int k=0;k<dim;++k) s += std::conj(F(k,i)) * Gz[static_cast<size_t>(k)];
                    b[static_cast<size_t>(i)] = s;
                }

                std::vector<complexe> Hz(dim, complexe(0.0));
                FullMatrix<complexe> A_copy = A; // FullMatrix::solve modifies a copy internally, but we keep simple
                A_copy.solve(Hz, b);
                const double nrm = vector_norm2(Hz);
                const double val = std::log(1.0 / (nrm + 1e-300));
                out << z1 << " " << z2 << " " << val << "\n";
            }
        }
    }

    int modes_max() const { return N; }

private:
    int N;
    int nm;
    double k0;
    double h;
    double L;
    Mode mode;
    std::vector<complexe> betas;
    FullMatrix<complexe> F;

    static double vector_norm2(const std::vector<complexe>& v){
        double s = 0.0;
        for(const auto& z : v) s += std::norm(z);
        return std::sqrt(s);
    }

    // Right-hand side G(z) from the statement (section 2.2).
    std::vector<complexe> build_rhs_G(double z1, double z2) const {
        const complexe ii = I_unit();

        if(mode == Mode::BackScattering){
            std::vector<complexe> G(static_cast<size_t>(nm), complexe(0.0));
            for(int m=0;m<nm;++m){
                const complexe b = betas[static_cast<size_t>(m)];
                const double cm = c_n(m, z2, h);
                G[static_cast<size_t>(m)] = std::exp(ii*b*(L + z1)) / (ii*b) * cm;
            }
            return G;
        }

        std::vector<complexe> G(static_cast<size_t>(2*nm), complexe(0.0));
        for(int m=0;m<nm;++m){
            const complexe b = betas[static_cast<size_t>(m)];
            const double cm = c_n(m, z2, h);
            G[static_cast<size_t>(m)]      = std::exp(ii*b*(L + z1)) / (ii*b) * cm;
            G[static_cast<size_t>(nm+m)]   = std::exp(ii*b*(L - z1)) / (ii*b) * cm;
        }
        return G;
    }
};

} // namespace usim

#endif // LINEAR_SAMPLING_HPP
