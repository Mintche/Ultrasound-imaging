#include "linear_sampling.hpp"

namespace LinearSampling {
    
#define EPSILON_LSM 1e-6 // Paramètre de régularisation Tikhonov

void compute_boundary_u_s(const MeshP2& mesh, int n_mode,int tag_left,int tag_right, 
                        double direction, [[maybe_unused]] double L, double k0, double h, 
                        std::vector<complexe>& u_s, const std::vector<complexe>& u_n) {

    for (const auto& nodes : mesh.nodes) {
        if (nodes.ref == tag_left || nodes.ref == tag_right){
            double c_n_val = Fem::evaluate_c_1d(nodes.y,h,n_mode);
            complexe b_n = Fem::compute_beta(k0,h,n_mode);
            u_s[nodes.id] = u_n[nodes.id] - c_n_val * exp(direction * complexe(0., 1.) * b_n * nodes.x);
        }
    }
}

void compute_projection(const std::vector<complexe>& u_s, const FullMatrix<complexe>& E_plus,
                        const FullMatrix<complexe>& E_minus, std::vector<complexe>& U_proj_plus, 
                        std::vector<complexe>& U_proj_minus) {

    U_proj_plus = E_plus.transpose() * u_s; // Projection sur les modes de droite
    U_proj_minus = E_minus.transpose() * u_s; // Projection sur les modes de gauche

}

FullMatrix<complexe> compute_F(const FullMatrix<complexe>& S_LL, 
                                const FullMatrix<complexe>& S_RL, 
                                const FullMatrix<complexe>& S_LR, 
                                const FullMatrix<complexe>& S_RR, 
                                int N_modes, double k0, double h, double L) {

    FullMatrix<complexe> F(2*N_modes, 2*N_modes);

    for (int n = 0; n < N_modes; ++n) {
        
        // Calcul du terme de phase et de normalisation (Eq. 231)
        complexe beta_n = Fem::compute_beta(k0, h, n);
        complexe i_beta = complexe(0., 1.) * beta_n;
        
        // Constante = e^(i * beta * L) / (i * beta)
        complexe constante = exp(i_beta * L) / i_beta;

        for (int m = 0; m < N_modes; ++m) {
            // Remplissage selon Eq (7) et définitions (231-235)
            
            // Bloc F++ (Haut-Gauche) : Utilise (Un-)+ => Source Droite, Mesure Droite
            F(m, n) = constante * S_RR(n, m); 
            
            // Bloc F+- (Haut-Droite) : Utilise (Un+)+ => Source Gauche, Mesure Droite
            F(m, n + N_modes) = constante * S_RL(n, m);
            
            // Bloc F-+ (Bas-Gauche)  : Utilise (Un-)- => Source Droite, Mesure Gauche
            F(m + N_modes, n) = constante * S_LR(n, m);
            
            // Bloc F-- (Bas-Droite)  : Utilise (Un+)- => Source Gauche, Mesure Gauche
            F(m + N_modes, n + N_modes) = constante * S_LL(n, m);
        }
    }
    return F;
}

FullMatrix<complexe> compute_F_pp(const FullMatrix<complexe>& S_m_p, int N_modes, double k0, double h, double L) {

    FullMatrix<complexe> F_pp(2*N_modes,2*N_modes);

    for (int n = 0; n < N_modes; ++n) {

        complexe exposant = complexe(0.,1.) * Fem::compute_beta(k0, h, n) * L; // exp(i*beta_n*L)
        complexe constante = exp(exposant) / (complexe(0.,1.) * Fem::compute_beta(k0, h, n)); // exp(i*beta_n*L) / (i*beta_n)

        for (int m = 0; m < N_modes; ++m) {
            F_pp(m, n) = constante * S_m_p(n,m); //F++ U-+
        }
    }
    return F_pp;
}

std::vector<complexe> assemble_Gz([[maybe_unused]] const MeshP2& mesh, int n_modes, double z1, double z2,
                                    double x_min, double x_max, double k0, double h){
    
    std::vector<complexe> Gz(2*n_modes);
    
    for (int m = 0; m < n_modes;m++){
        
        complexe denominateur = complexe(0.,1.) * Fem::compute_beta(k0, h, m);
        complexe c_val = Fem::evaluate_c_1d(z2, h, m); // Attention: suppose y dans [0, h]
        
        // Phase relative au bord gauche (x_min) et droit (x_max)
        Gz[m] = exp(complexe(0.,1.) * Fem::compute_beta(k0, h, m) * (x_max - z1)) * c_val / denominateur;
        Gz[m + n_modes] = exp(complexe(0.,1.) * Fem::compute_beta(k0, h, m) * (z1 - x_min)) * c_val / denominateur;
    }
    return Gz;
}

void add_gaussian_noise(std::vector<complexe>& data, double noise_level) {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> distribution(0.0, noise_level);

    for (auto& val : data) {
        double noise_real = distribution(generator);
        double noise_imag = distribution(generator);
        val += complexe(noise_real, noise_imag);
    }
}

void compute_lsm_single_freq(const MeshP2& mesh, double k0, double kd, double noise_level,
                             int grid_nx, int grid_ny, 
                             double x_scan_min, double x_scan_max, double y_scan_min, double y_scan_max,
                             int tag_left, int tag_right, std::vector<double>& indicators) {

    double h = mesh.Ly;
    double L = mesh.Lx / 2.0;
    double x_source_gauche = mesh.xmin;
    double x_source_droite = mesh.xmax;

    // Calcul de N_MODES (nombre de modes de guide d'ondes)
    int N_MODES = floor(h * std::min(k0,kd) / M_PI) + 5; 

    printf("  -> Calcul k0=%.2f, kd=%.2f | N_modes=%d\n", k0, kd, N_MODES);

    // Calcul du profil et matrices FEM
    std::vector<size_t> profile = Fem::compute_profile_enhanced(mesh, {tag_left, tag_right});
    ProfileMatrix<complexe> K(profile);

    Fem::A_matrix(mesh, K, 1.0);                 // Matrice de Raideur
    Fem::B_matrix(mesh, K, k0, kd, -1.0);        // Matrice de Masse
    
    // Matrices de projection et diagonale D
    FullMatrix<complexe> E_minus = Fem::compute_E(mesh, N_MODES, tag_left, k0);
    FullMatrix<complexe> E_plus  = Fem::compute_E(mesh, N_MODES, tag_right, k0);
    FullMatrix<complexe> D(N_MODES, N_MODES);
    Fem::compute_D(D, N_MODES, h, k0);

    // Conditions DtN
    Fem::T_matrix(K, E_minus, D, h, tag_left, -1.0);
    Fem::T_matrix(K, E_plus, D, h, tag_right, -1.0);

    // Factorisation
    K.factorize(); 

    // Simulation et construction des matrices S
    FullMatrix<complexe> S_LL(N_MODES, N_MODES);
    FullMatrix<complexe> S_RL(N_MODES, N_MODES);
    FullMatrix<complexe> S_LR(N_MODES, N_MODES);
    FullMatrix<complexe> S_RR(N_MODES, N_MODES);

    std::vector<complexe> U(mesh.ndof());
    std::vector<complexe> u_s(mesh.ndof());
    std::vector<complexe> proj_plus(N_MODES), proj_minus(N_MODES);

    // Cas 1 : Source Gauche
    for (int n = 0; n < N_MODES; ++n) {
        std::vector<complexe> G = Fem::assemble_source_vector(mesh, E_minus, n, k0, x_source_gauche, 1.0);
        K.solve(U, G); 
        if (noise_level > 0) LinearSampling::add_gaussian_noise(U, noise_level);
        LinearSampling::compute_boundary_u_s(mesh, n, tag_left, tag_right, 1.0, x_source_gauche, k0, h, u_s, U);
        LinearSampling::compute_projection(u_s, E_plus, E_minus, proj_plus, proj_minus);
        for(int m=0; m<N_MODES; ++m) {
            S_RL(n, m) = proj_plus[m];
            S_LL(n, m) = proj_minus[m];
        }
    }

    // Cas 2 : Source Droite
    for (int n = 0; n < N_MODES; ++n) {
        std::vector<complexe> G = Fem::assemble_source_vector(mesh, E_plus, n, k0, x_source_droite, -1.0); 
        K.solve(U, G);
        if (noise_level > 0) LinearSampling::add_gaussian_noise(U, noise_level);
        LinearSampling::compute_boundary_u_s(mesh, n, tag_left, tag_right, -1.0, x_source_droite, k0, h, u_s, U);
        LinearSampling::compute_projection(u_s, E_plus, E_minus, proj_plus, proj_minus);
        for(int m=0; m<N_MODES; ++m) {
            S_RR(n, m) = proj_plus[m];
            S_LR(n, m) = proj_minus[m];
        }
    }

    // Construction LSM
    FullMatrix<complexe> F = LinearSampling::compute_F(S_LL, S_RL, S_LR, S_RR, N_MODES, k0, h, L);
    FullMatrix<complexe> F_adj = F.adjoint();
    FullMatrix<complexe> I(2*N_MODES,2*N_MODES);
    for (int i = 0; i < 2*N_MODES;i++) I(i,i) = 1.0;
    FullMatrix<complexe> M = F_adj * F + std::complex(EPSILON_LSM, 0.0)*I;
    M.factorize();

    // Calcul sur la grille
    indicators.resize(grid_nx * grid_ny);
    for(int i = 0; i < grid_nx; ++i) {
        for(int j = 0; j < grid_ny; ++j) {
            std::vector<complexe> H_z(2*N_MODES);
            double z1 = x_scan_min + i * (x_scan_max - x_scan_min) / (grid_nx - 1);
            double z2 = y_scan_min + j * (y_scan_max - y_scan_min) / (grid_ny - 1);

            std::vector<complexe> G_z = LinearSampling::assemble_Gz(mesh, N_MODES, z1, z2 - mesh.ymin, mesh.xmin, mesh.xmax, k0, h);
            std::vector<complexe> F_adj_G = F_adj * G_z;
            M.solve(H_z, F_adj_G);

            double norm_Hz = norm(H_z); 
            indicators[i * grid_ny + j] = 1.0 / (norm_Hz + 1e-15);
        }
    }
}

void compute_lsm_average(const MeshP2& mesh, int n_freq, double base_k0, double contrast_ratio,
                         double noise_percentage, int grid_nx, int grid_ny,
                         int tag_left, int tag_right, const std::string& output_filename) {

    double h = mesh.Ly;
    double noise_level = (noise_percentage <= 0) ? -1 : noise_percentage * (1.0/sqrt(2.0*h));
    
    // Marges de sécurité pour le scan
    double x_scan_min = mesh.xmin + 0.05; double x_scan_max = mesh.xmax - 0.05;
    double y_scan_min = mesh.ymin + 0.05; double y_scan_max = mesh.ymax - 0.05;

    std::vector<double> avg_indicators(grid_nx * grid_ny, 0.0);
    std::vector<double> current_indicators;

    for(int i = 0; i < n_freq; ++i) {
        double k0 = base_k0 + i * 1.0; // Incrément de 1.0 par pas de fréquence
        double kd = k0 * contrast_ratio;
        
        compute_lsm_single_freq(mesh, k0, kd, noise_level, grid_nx, grid_ny, 
                                x_scan_min, x_scan_max, y_scan_min, y_scan_max, 
                                tag_left, tag_right, current_indicators);

        double max_val = normesup(current_indicators);
        for(size_t k=0; k<avg_indicators.size(); ++k) {
            avg_indicators[k] += (current_indicators[k] / max_val);
        }
    }

    std::ofstream file(output_filename);
    file << grid_nx << " " << grid_ny << "\n";
    for(int i = 0; i < grid_nx; ++i) {
        double z1 = x_scan_min + i * (x_scan_max - x_scan_min) / (grid_nx - 1);
        for(int j = 0; j < grid_ny; ++j) {
            double z2 = y_scan_min + j * (y_scan_max - y_scan_min) / (grid_ny - 1);
            file << z1 << " " << z2 << " " << (avg_indicators[i * grid_ny + j] / n_freq) << "\n";
        }
    }
    file.close();
    printf("Resultats ecrits dans '%s'.\n", output_filename.c_str());
}

void compute_lsm_physical(const MeshP2& mesh, PhysicalParameters phys_params, double contrast_ratio,
                          double noise_percentage, int grid_nx, int grid_ny,
                          int tag_left, int tag_right, const std::string& output_filename) {

    // 1. Analyse du maillage
    double h_max = Fem::get_max_edge_length(mesh);
    printf("--- Verification du maillage ---\n");
    printf("Pas maximum du maillage (h_max) : %.4f m\n", h_max);

    // 2. Vérification pour la fréquence maximale
    double freq_max = phys_params.freq_start + (phys_params.n_freq - 1) * phys_params.freq_step;
    double lambda_min = phys_params.c0 / freq_max;
    
    // Critère P2 : On recommande h < lambda / 3 (environ 6 noeuds par longueur d'onde)
    // Le noeud milieu aide, donc h correspond à 2 intervalles nodaux.
    double points_per_wavelength = lambda_min / (h_max / 2.0);

    printf("Frequence max : %.2f Hz => Longueur d'onde min : %.4f m\n", freq_max, lambda_min);
    printf("Points par longueur d'onde (approx) : %.2f\n", points_per_wavelength);

    if (points_per_wavelength < 6.0) {
        printf("\n[ATTENTION] Le maillage est trop grossier pour la frequence demandee !\n");
        printf("Conseil : Reduire h_max en dessous de %.4f m ou baisser la frequence.\n\n", lambda_min / 3.0);
        throw std::runtime_error("Maillage inadapte pour la frequence.");
    } else {
        printf("[OK] Resolution du maillage suffisante.\n\n");
    }

    // 3. Conversion et appel de la fonction générique
    // k = 2*pi*f / c0

    double h = mesh.Ly;
    double noise_level = (noise_percentage <= 0) ? -1 : noise_percentage * (1.0/sqrt(2.0*h));
    
    double x_scan_min = mesh.xmin + 0.05; double x_scan_max = mesh.xmax - 0.05;
    double y_scan_min = mesh.ymin + 0.05; double y_scan_max = mesh.ymax - 0.05;

    std::vector<double> avg_indicators(grid_nx * grid_ny, 0.0);
    std::vector<double> current_indicators;

    for(int i = 0; i < phys_params.n_freq; ++i) {
        double f = phys_params.freq_start + i * phys_params.freq_step;
        double k0 = 2.0 * M_PI * f / phys_params.c0;
        double kd = k0 * contrast_ratio;
        
        printf("Simulation Freq=%.1f Hz (k0=%.2f) ...\n", f, k0);

        compute_lsm_single_freq(mesh, k0, kd, noise_level, grid_nx, grid_ny, 
                                x_scan_min, x_scan_max, y_scan_min, y_scan_max, 
                                tag_left, tag_right, current_indicators);

        double max_val = normesup(current_indicators);
        for(size_t k=0; k<avg_indicators.size(); ++k) {
            avg_indicators[k] += (current_indicators[k] / max_val);
        }
    }

    // Ecriture fichier
    std::ofstream file(output_filename);
    file << grid_nx << " " << grid_ny << "\n";
    for(int i = 0; i < grid_nx; ++i) {
        double z1 = x_scan_min + i * (x_scan_max - x_scan_min) / (grid_nx - 1);
        for(int j = 0; j < grid_ny; ++j) {
            double z2 = y_scan_min + j * (y_scan_max - y_scan_min) / (grid_ny - 1);
            file << z1 << " " << z2 << " " << (avg_indicators[i * grid_ny + j] / phys_params.n_freq) << "\n";
        }
    }
    file.close();
    printf("Resultats (Physique) ecrits dans '%s'.\n", output_filename.c_str());
}

} // namespace LinearSampling