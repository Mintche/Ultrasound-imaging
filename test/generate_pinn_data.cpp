#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "fem.hpp"
#include "p2_interpolation.hpp"

using namespace std;
using namespace usim;

namespace {

struct Configuration {
    filesystem::path mesh_file;
    filesystem::path output_dir = "pinn_data";
    string dataset_name;
    double c0 = 340.0;
    double defect_speed_ratio = 1.0;
    vector<double> frequencies;
    vector<int> modes{0};
    int grid_nx = 201;
    int grid_ny = 61;
    bool export_grid = true;
    int defect_tag = 2;
    int left_tag = 11;
    int right_tag = 12;
};

struct GridSample {
    double x = 0.0;
    double y = 0.0;
    P2PointLocation location;
};

struct ModeFiles {
    ofstream left;
    ofstream right;
    ofstream field;
    ofstream grid;
};

void print_usage(const char* executable) {
    cout
        << "Usage:\n  " << executable << " --mesh mesh.msh --dataset NAME"
        << " --frequencies F1,F2 --modes M1,M2 [options]\n\n"
        << "Options:\n"
        << "  --output-dir DIR       Dossier de sortie (defaut: pinn_data)\n"
        << "  --c0 VALUE             Celerite du milieu sain (defaut: 340)\n"
        << "  --speed-ratio VALUE    c_defaut / c0 (defaut: 1)\n"
        << "  --grid NXxNY           Grille reguliere interpolee (defaut: 201x61)\n"
        << "  --no-grid              Ne pas exporter la grille interpolee\n"
        << "  --defect-tag TAG       Tag physique du defaut (defaut: 2)\n"
        << "  --left-tag TAG         Tag du port gauche (defaut: 11)\n"
        << "  --right-tag TAG        Tag du port droit (defaut: 12)\n"
        << "  --help                  Afficher cette aide\n\n"
        << "Exemple:\n  " << executable
        << " --mesh data/test_us_barrehalf_centree.msh"
        << " --output-dir pinn_data --dataset barrehalf_contrast20percent"
        << " --c0 340 --speed-ratio 0.8 --frequencies 500,600,700"
        << " --modes 0,1 --grid 201x61\n\n"
        << "L'ancien format positionnel reste accepte:\n  " << executable
        << " mesh.msh mode c0 speed_ratio f1 [f2 ...]\n";
}

template <typename T, typename Converter>
vector<T> parse_list(const string& text, Converter converter, const string& option) {
    vector<T> values;
    string token;
    stringstream stream(text);
    while (getline(stream, token, ',')) {
        if (token.empty()) {
            throw invalid_argument("Valeur vide dans " + option);
        }
        values.push_back(converter(token));
    }
    if (values.empty()) {
        throw invalid_argument("La liste " + option + " est vide");
    }
    return values;
}

string sanitize_name(string name) {
    for (char& character : name) {
        const bool valid = (character >= 'a' && character <= 'z') ||
                           (character >= 'A' && character <= 'Z') ||
                           (character >= '0' && character <= '9') ||
                           character == '-' || character == '_';
        if (!valid) character = '_';
    }
    if (name.empty()) name = "dataset";
    return name;
}

Configuration parse_arguments(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        throw invalid_argument("Arguments manquants");
    }

    Configuration config;
    const string first_argument = argv[1];

    // Compatibilite avec l'interface initiale du projet.
    if (!first_argument.empty() && first_argument[0] != '-') {
        if (argc < 6) {
            print_usage(argv[0]);
            throw invalid_argument("Le format positionnel attend au moins 5 arguments");
        }
        config.mesh_file = argv[1];
        config.modes = {stoi(argv[2])};
        config.c0 = stod(argv[3]);
        config.defect_speed_ratio = stod(argv[4]);
        for (int i = 5; i < argc; ++i) config.frequencies.push_back(stod(argv[i]));
        config.dataset_name = config.mesh_file.stem().string();
        config.output_dir = ".";
        cerr << "Attention: interface positionnelle historique utilisee. "
             << "Utilisez --help pour la nouvelle interface.\n";
        return config;
    }

    auto value_after = [&](int& index, const string& option) -> string {
        if (index + 1 >= argc) throw invalid_argument("Valeur manquante apres " + option);
        return argv[++index];
    };

    for (int i = 1; i < argc; ++i) {
        const string option = argv[i];
        if (option == "--help" || option == "-h") {
            print_usage(argv[0]);
            exit(0);
        } else if (option == "--mesh") {
            config.mesh_file = value_after(i, option);
        } else if (option == "--output-dir") {
            config.output_dir = value_after(i, option);
        } else if (option == "--dataset") {
            config.dataset_name = value_after(i, option);
        } else if (option == "--c0") {
            config.c0 = stod(value_after(i, option));
        } else if (option == "--speed-ratio") {
            config.defect_speed_ratio = stod(value_after(i, option));
        } else if (option == "--frequencies") {
            config.frequencies = parse_list<double>(
                value_after(i, option), [](const string& value) { return stod(value); }, option);
        } else if (option == "--modes") {
            config.modes = parse_list<int>(
                value_after(i, option), [](const string& value) { return stoi(value); }, option);
        } else if (option == "--grid") {
            const string dimensions = value_after(i, option);
            const size_t separator = dimensions.find_first_of("xX");
            if (separator == string::npos) throw invalid_argument("--grid attend NXxNY");
            config.grid_nx = stoi(dimensions.substr(0, separator));
            config.grid_ny = stoi(dimensions.substr(separator + 1));
            config.export_grid = true;
        } else if (option == "--no-grid") {
            config.export_grid = false;
        } else if (option == "--defect-tag") {
            config.defect_tag = stoi(value_after(i, option));
        } else if (option == "--left-tag") {
            config.left_tag = stoi(value_after(i, option));
        } else if (option == "--right-tag") {
            config.right_tag = stoi(value_after(i, option));
        } else {
            throw invalid_argument("Option inconnue: " + option);
        }
    }

    if (config.mesh_file.empty()) throw invalid_argument("--mesh est obligatoire");
    if (config.dataset_name.empty()) config.dataset_name = config.mesh_file.stem().string();
    config.dataset_name = sanitize_name(config.dataset_name);
    if (config.frequencies.empty()) throw invalid_argument("--frequencies est obligatoire");
    if (!(config.c0 > 0.0)) throw invalid_argument("c0 doit etre strictement positif");
    if (!(config.defect_speed_ratio > 0.0)) {
        throw invalid_argument("speed-ratio doit etre strictement positif");
    }
    if (config.grid_nx < 2 || config.grid_ny < 2) {
        throw invalid_argument("La grille doit contenir au moins 2x2 points");
    }
    for (double frequency : config.frequencies) {
        if (!(frequency > 0.0)) throw invalid_argument("Toutes les frequences doivent etre positives");
    }
    for (int mode : config.modes) {
        if (mode < 0) throw invalid_argument("Les indices de mode doivent etre positifs ou nuls");
    }
    sort(config.frequencies.begin(), config.frequencies.end());
    config.frequencies.erase(unique(config.frequencies.begin(), config.frequencies.end()),
                             config.frequencies.end());
    sort(config.modes.begin(), config.modes.end());
    config.modes.erase(unique(config.modes.begin(), config.modes.end()), config.modes.end());
    return config;
}

ofstream open_output(const filesystem::path& path) {
    ofstream output(path);
    if (!output) throw runtime_error("Impossible de creer " + path.string());
    output << setprecision(17);
    return output;
}

vector<int> boundary_nodes(const MeshP2& mesh, int boundary_tag) {
    set<int> unique_nodes;
    for (const auto& triangle : mesh.triangles) {
        for (int edge = 0; edge < 3; ++edge) {
            if (triangle.edge_ref[edge] != boundary_tag) continue;
            unique_nodes.insert(triangle.node_ids[edge]);
            unique_nodes.insert(triangle.node_ids[(edge + 1) % 3]);
            unique_nodes.insert(triangle.node_ids[edge + 3]);
        }
    }
    vector<int> nodes(unique_nodes.begin(), unique_nodes.end());
    sort(nodes.begin(), nodes.end(), [&](int first, int second) {
        const auto& a = mesh.nodes[first];
        const auto& b = mesh.nodes[second];
        if (a.y != b.y) return a.y < b.y;
        return a.x < b.x;
    });
    return nodes;
}

vector<GridSample> prepare_grid(const MeshP2& mesh, int nx, int ny) {
    const size_t number_of_samples = static_cast<size_t>(nx) * static_cast<size_t>(ny);
    vector<vector<size_t>> candidate_triangles(number_of_samples);

    // Spatial index specialized for the requested regular grid. Each triangle is
    // registered only for grid points contained in its bounding box.
    for (size_t triangle_id = 0; triangle_id < mesh.triangles.size(); ++triangle_id) {
        const auto& triangle = mesh.triangles[triangle_id];
        const auto& p0 = mesh.nodes[triangle.node_ids[0]];
        const auto& p1 = mesh.nodes[triangle.node_ids[1]];
        const auto& p2 = mesh.nodes[triangle.node_ids[2]];
        const double min_x = min({p0.x, p1.x, p2.x});
        const double max_x = max({p0.x, p1.x, p2.x});
        const double min_y = min({p0.y, p1.y, p2.y});
        const double max_y = max({p0.y, p1.y, p2.y});

        const double scaled_x_min = (min_x - mesh.xmin) * (nx - 1) / mesh.Lx;
        const double scaled_x_max = (max_x - mesh.xmin) * (nx - 1) / mesh.Lx;
        const double scaled_y_min = (min_y - mesh.ymin) * (ny - 1) / mesh.Ly;
        const double scaled_y_max = (max_y - mesh.ymin) * (ny - 1) / mesh.Ly;
        const int ix_min = max(0, static_cast<int>(ceil(scaled_x_min - 1e-10)));
        const int ix_max = min(nx - 1, static_cast<int>(floor(scaled_x_max + 1e-10)));
        const int iy_min = max(0, static_cast<int>(ceil(scaled_y_min - 1e-10)));
        const int iy_max = min(ny - 1, static_cast<int>(floor(scaled_y_max + 1e-10)));

        for (int iy = iy_min; iy <= iy_max; ++iy) {
            for (int ix = ix_min; ix <= ix_max; ++ix) {
                candidate_triangles[static_cast<size_t>(iy) * nx + ix].push_back(triangle_id);
            }
        }
    }

    vector<GridSample> samples;
    samples.reserve(number_of_samples);
    for (int iy = 0; iy < ny; ++iy) {
        const double y = mesh.ymin + mesh.Ly * static_cast<double>(iy) / (ny - 1);
        for (int ix = 0; ix < nx; ++ix) {
            const double x = mesh.xmin + mesh.Lx * static_cast<double>(ix) / (nx - 1);
            P2PointLocation location;
            for (size_t triangle_id : candidate_triangles[static_cast<size_t>(iy) * nx + ix]) {
                location = locate_point_in_triangle_p2(mesh, triangle_id, x, y);
                if (location.found()) break;
            }
            if (!location.found()) {
                throw runtime_error("Point de grille hors maillage: (" + to_string(x) + ", " +
                                    to_string(y) + ")");
            }
            samples.push_back({x, y, location});
        }
    }
    return samples;
}

void export_mesh(const Configuration& config, const MeshP2& mesh) {
    auto nodes = open_output(config.output_dir /
                             ("fem_mesh_nodes_" + config.dataset_name + ".csv"));
    nodes << "node_id,x,y,x_norm,y_norm,boundary_tag\n";
    for (const auto& node : mesh.nodes) {
        const double x_norm = 2.0 * (node.x - mesh.xmin) / mesh.Lx - 1.0;
        const double y_norm = 2.0 * (node.y - mesh.ymin) / mesh.Ly - 1.0;
        nodes << node.id << ',' << node.x << ',' << node.y << ',' << x_norm << ',' << y_norm
              << ',' << node.ref << '\n';
    }

    auto elements = open_output(config.output_dir /
                                ("fem_mesh_elements_" + config.dataset_name + ".csv"));
    elements << "element_id,node0,node1,node2,node3,node4,node5,region_tag,is_defect\n";
    for (size_t element_id = 0; element_id < mesh.triangles.size(); ++element_id) {
        const auto& triangle = mesh.triangles[element_id];
        elements << element_id;
        for (int node_id : triangle.node_ids) elements << ',' << node_id;
        elements << ',' << triangle.ref << ',' << static_cast<int>(triangle.is_defect) << '\n';
    }
}

void export_evaluation_grid(const Configuration& config, const MeshP2& mesh,
                            const vector<GridSample>& samples) {
    auto output = open_output(config.output_dir /
                              ("fem_evaluation_grid_" + config.dataset_name + ".csv"));
    output << "grid_i,x,y,x_norm,y_norm,element_id,xi,eta,"
              "node0,node1,node2,node3,node4,node5,"
              "phi0,phi1,phi2,phi3,phi4,phi5,region_tag,is_defect,c\n";

    for (size_t grid_i = 0; grid_i < samples.size(); ++grid_i) {
        const auto& sample = samples[grid_i];
        const auto& triangle = mesh.triangles[sample.location.triangle_id];
        const auto phi = p2_shape_functions(sample.location.xi, sample.location.eta);
        const double x_norm = 2.0 * (sample.x - mesh.xmin) / mesh.Lx - 1.0;
        const double y_norm = 2.0 * (sample.y - mesh.ymin) / mesh.Ly - 1.0;
        const double sound_speed =
            config.c0 * (triangle.is_defect ? config.defect_speed_ratio : 1.0);

        output << grid_i << ',' << sample.x << ',' << sample.y << ',' << x_norm << ','
               << y_norm << ',' << sample.location.triangle_id << ','
               << sample.location.xi << ',' << sample.location.eta;
        for (int node_id : triangle.node_ids) output << ',' << node_id;
        for (double weight : phi) output << ',' << weight;
        output << ',' << triangle.ref << ',' << static_cast<int>(triangle.is_defect) << ','
               << sound_speed << '\n';
    }
}

ModeFiles open_mode_files(const Configuration& config, int mode) {
    const string suffix = config.dataset_name + "_mode" + to_string(mode) + ".csv";
    ModeFiles files{
        open_output(config.output_dir / ("pinn_boundary_left_" + suffix)),
        open_output(config.output_dir / ("pinn_boundary_right_" + suffix)),
        open_output(config.output_dir / ("fem_field_" + suffix)),
        ofstream{}
    };
    files.left << "f,k0,x,y,Re_U,Im_U\n";
    files.right << "f,k0,x,y,Re_U,Im_U\n";
    files.field << "f,k0,mode,node_id,x,y,x_norm,y_norm,Re_U,Im_U,abs_U\n";
    if (config.export_grid) {
        files.grid = open_output(config.output_dir / ("fem_grid_" + suffix));
        files.grid << "f,k0,mode,grid_i,x,y,x_norm,y_norm,element_id,region_tag,is_defect,c,Re_U,Im_U,abs_U\n";
    }
    return files;
}

void export_boundary(ofstream& output, const MeshP2& mesh, const vector<int>& node_ids,
                     const vector<complexe>& field, double frequency, double k0) {
    for (int node_id : node_ids) {
        const auto& node = mesh.nodes[node_id];
        output << frequency << ',' << k0 << ',' << node.x << ',' << node.y << ','
               << real(field[node_id]) << ',' << imag(field[node_id]) << '\n';
    }
}

void export_nodal_field(ofstream& output, const MeshP2& mesh, const vector<complexe>& field,
                        double frequency, double k0, int mode) {
    for (const auto& node : mesh.nodes) {
        const double x_norm = 2.0 * (node.x - mesh.xmin) / mesh.Lx - 1.0;
        const double y_norm = 2.0 * (node.y - mesh.ymin) / mesh.Ly - 1.0;
        output << frequency << ',' << k0 << ',' << mode << ',' << node.id << ',' << node.x
               << ',' << node.y << ',' << x_norm << ',' << y_norm << ','
               << real(field[node.id]) << ',' << imag(field[node.id]) << ','
               << abs(field[node.id]) << '\n';
    }
}

void export_grid_field(ofstream& output, const MeshP2& mesh,
                       const vector<GridSample>& samples, const vector<complexe>& field,
                       double frequency, double k0, int mode, double c0,
                       double defect_speed_ratio) {
    for (size_t grid_i = 0; grid_i < samples.size(); ++grid_i) {
        const auto& sample = samples[grid_i];
        const auto& triangle = mesh.triangles[sample.location.triangle_id];
        const complexe value = interpolate_p2(mesh, field, sample.location);
        const double x_norm = 2.0 * (sample.x - mesh.xmin) / mesh.Lx - 1.0;
        const double y_norm = 2.0 * (sample.y - mesh.ymin) / mesh.Ly - 1.0;
        const double sound_speed = c0 * (triangle.is_defect ? defect_speed_ratio : 1.0);
        output << frequency << ',' << k0 << ',' << mode << ',' << grid_i << ',' << sample.x
               << ',' << sample.y << ',' << x_norm << ',' << y_norm << ','
               << sample.location.triangle_id << ',' << triangle.ref << ','
               << static_cast<int>(triangle.is_defect) << ',' << sound_speed << ','
               << real(value) << ',' << imag(value) << ',' << abs(value) << '\n';
    }
}

void export_metadata(const Configuration& config, const MeshP2& mesh) {
    auto metadata = open_output(config.output_dir /
                                ("fem_metadata_" + config.dataset_name + ".txt"));
    metadata << "schema_version=2\n"
             << "mesh=" << config.mesh_file.string() << '\n'
             << "dataset=" << config.dataset_name << '\n'
             << "element=P2_triangle\n"
             << "node_index_base=0\n"
             << "c0=" << config.c0 << '\n'
             << "defect_speed_ratio=" << config.defect_speed_ratio << '\n'
             << "defect_tag=" << config.defect_tag << '\n'
             << "left_tag=" << config.left_tag << '\n'
             << "right_tag=" << config.right_tag << '\n';
    metadata << "frequencies=";
    for (size_t i = 0; i < config.frequencies.size(); ++i) {
        if (i > 0) metadata << ',';
        metadata << config.frequencies[i];
    }
    metadata << "\nmodes=";
    for (size_t i = 0; i < config.modes.size(); ++i) {
        if (i > 0) metadata << ',';
        metadata << config.modes[i];
    }
    metadata << "\nxmin=" << mesh.xmin << '\n'
             << "xmax=" << mesh.xmax << '\n'
             << "ymin=" << mesh.ymin << '\n'
             << "ymax=" << mesh.ymax << '\n'
             << "ndof=" << mesh.ndof() << '\n'
             << "n_elements=" << mesh.triangles.size() << '\n'
             << "grid_nx=" << (config.export_grid ? config.grid_nx : 0) << '\n'
             << "grid_ny=" << (config.export_grid ? config.grid_ny : 0) << '\n'
             << "grid_order=x_fastest\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Configuration config = parse_arguments(argc, argv);
        filesystem::create_directories(config.output_dir);

        MeshP2 mesh;
        mesh.read_msh_v2_ascii(config.mesh_file.string(), {config.defect_tag});
        Fem::reorder_mesh_rcm(mesh);

        if (!(mesh.Lx > 0.0 && mesh.Ly > 0.0)) {
            throw runtime_error("Le maillage doit avoir des dimensions strictement positives");
        }

        const vector<int> left_nodes = boundary_nodes(mesh, config.left_tag);
        const vector<int> right_nodes = boundary_nodes(mesh, config.right_tag);
        if (left_nodes.empty() || right_nodes.empty()) {
            throw runtime_error("Les tags de port ne correspondent a aucune arete du maillage");
        }

        cout << "Maillage: " << mesh.ndof() << " ddl P2, " << mesh.triangles.size()
             << " triangles, ports " << left_nodes.size() << "/" << right_nodes.size()
             << " noeuds.\n";

        export_mesh(config, mesh);
        export_metadata(config, mesh);

        vector<GridSample> grid_samples;
        if (config.export_grid) {
            cout << "Localisation des " << config.grid_nx * config.grid_ny
                 << " points de la grille...\n";
            grid_samples = prepare_grid(mesh, config.grid_nx, config.grid_ny);
            export_evaluation_grid(config, mesh, grid_samples);
        }

        vector<ModeFiles> mode_files;
        mode_files.reserve(config.modes.size());
        for (int mode : config.modes) mode_files.push_back(open_mode_files(config, mode));

        const vector<size_t> profile =
            Fem::compute_profile_enhanced(mesh, {config.left_tag, config.right_tag});
        const int highest_requested_mode = *max_element(config.modes.begin(), config.modes.end());

        for (double frequency : config.frequencies) {
            const double k0 = 2.0 * M_PI * frequency / config.c0;
            const double kd = k0 / config.defect_speed_ratio;
            const int number_of_dtn_modes =
                max(static_cast<int>(floor(mesh.Ly * k0 / M_PI)) + 5,
                    highest_requested_mode + 1);

            cout << "Resolution f=" << frequency << " Hz, k0=" << k0
                 << ", " << number_of_dtn_modes << " modes DtN...\n";

            ProfileMatrix<complexe> system(profile);
            Fem::A_matrix(mesh, system, 1.0);
            Fem::B_matrix(mesh, system, k0, kd, -1.0);

            FullMatrix<complexe> e_left =
                Fem::compute_E(mesh, number_of_dtn_modes, config.left_tag, k0);
            FullMatrix<complexe> e_right =
                Fem::compute_E(mesh, number_of_dtn_modes, config.right_tag, k0);
            FullMatrix<complexe> dtn(number_of_dtn_modes, number_of_dtn_modes);
            Fem::compute_D(dtn, number_of_dtn_modes, mesh.Ly, k0);
            Fem::T_matrix(system, e_left, dtn, mesh.Ly, config.left_tag, -1.0);
            Fem::T_matrix(system, e_right, dtn, mesh.Ly, config.right_tag, -1.0);
            system.factorize();

            for (size_t mode_index = 0; mode_index < config.modes.size(); ++mode_index) {
                const int mode = config.modes[mode_index];
                const complexe beta = Fem::compute_beta(k0, mesh.Ly, mode);
                if (abs(imag(beta)) > 1e-12) {
                    cerr << "Attention: le mode incident " << mode << " est evanescent a "
                         << frequency << " Hz.\n";
                }

                const vector<complexe> source =
                    Fem::assemble_source_vector(mesh, e_left, mode, k0, mesh.xmin, 1.0);
                vector<complexe> field(mesh.ndof());
                system.solve(field, source);

                auto& files = mode_files[mode_index];
                export_boundary(files.left, mesh, left_nodes, field, frequency, k0);
                export_boundary(files.right, mesh, right_nodes, field, frequency, k0);
                export_nodal_field(files.field, mesh, field, frequency, k0, mode);
                if (config.export_grid) {
                    export_grid_field(files.grid, mesh, grid_samples, field, frequency, k0,
                                      mode, config.c0, config.defect_speed_ratio);
                }
            }
        }

        cout << "Donnees FEM/PINN generees dans " << config.output_dir << "\n";
        return 0;
    } catch (const exception& error) {
        cerr << "Erreur: " << error.what() << '\n';
        return 1;
    }
}
