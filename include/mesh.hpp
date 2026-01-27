#ifndef MESH_HPP
#define MESH_HPP

#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "math.hpp" // for complexe (std::complex<double>)

/*
  Mesh reader + P2 enrichment (Student E2)

  - Reads a 2D triangular mesh produced by Gmsh (.msh).
  - Supports (reliably) Gmsh MSH v2 ASCII.
  - Gmsh often outputs only P1 triangles (3 vertices). For P2 FEM we must create
    the 3 mid-edge nodes and share them between adjacent triangles.
  - Also records boundary references on triangle edges from boundary line elements.

  Data model (per project statement):
  - Mesh: list of triangles.
  - Triangle: 6 nodes (3 vertices + 3 mid-edge points), a region/reference id,
    a list of global node numbers, and an edge-boundary marker (0 if not on boundary).
*/

namespace usim {

struct Point2D {
    double x = 0.0;
    double y = 0.0;
    int    id = -1;   // global index (0-based)
    int    ref = 0;   // optional node reference (0 if none)
};

struct TriangleP2 {
    // Local ordering compatible with the statement:
    // 1-2-3 are vertices, 4 is midpoint(1,2), 5 is midpoint(2,3), 6 is midpoint(1,3).
    std::array<int,6> node_ids{{-1,-1,-1,-1,-1,-1}};
    int ref = 0;                // triangle "physical" tag (e.g. defect zone id)
    bool is_defect = false;     // user can set from ref via set_defect_tags(...)
    std::array<int,3> edge_ref{{0,0,0}}; // for edges (1,2), (2,3), (3,1): boundary tag or 0
};

namespace detail {

// Trim helpers
inline std::string ltrim(std::string s){
    size_t i=0; while(i<s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    return s.substr(i);
}
inline std::string rtrim(std::string s){
    if(s.empty()) return s;
    size_t i=s.size();
    while(i>0 && std::isspace(static_cast<unsigned char>(s[i-1]))) --i;
    return s.substr(0,i);
}
inline std::string trim(std::string s){ return rtrim(ltrim(std::move(s))); }

inline bool starts_with(const std::string& s, const std::string& p){
    return s.size()>=p.size() && s.compare(0,p.size(),p)==0;
}

inline std::vector<long long> split_ll(const std::string& line){
    std::stringstream ss(line);
    std::vector<long long> out;
    long long v;
    while(ss>>v) out.push_back(v);
    return out;
}

// Edge key for unordered_map: store ordered (min,max)
struct EdgeKey {
    int a=0, b=0;
    EdgeKey() = default;
    EdgeKey(int i, int j){
        if(i<j){ a=i; b=j; }
        else { a=j; b=i; }
    }
    bool operator==(const EdgeKey& o) const { return a==o.a && b==o.b; }
};

struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& e) const noexcept {
        // simple hash combine for two ints
        return (static_cast<std::size_t>(e.a) << 32) ^ static_cast<std::size_t>(e.b);
    }
};

} // namespace detail

class MeshP2 {
public:
    // Public data for simplicity (project-style).
    std::vector<Point2D>    nodes;      // includes generated mid-edge nodes
    std::vector<TriangleP2> triangles;  // P2 triangles

    // Read a .msh (MSH v2 ASCII) and enrich it to P2.
    // defect_tags: if not empty, triangles with ref in defect_tags are marked as defect.
    void read_msh_v2_ascii(const std::string& filename,  const std::vector<int>& defect_tags = {}) {

        std::ifstream in(filename);
        if(!in) throw std::runtime_error("MeshP2: cannot open file: " + filename);

        // --- Parse $MeshFormat (optional in some exports but usually present) ---
        std::string line;
        double msh_version = 2.2;
        bool is_binary = false;

        while(std::getline(in,line)){
            line = detail::trim(line);
            if(line=="$MeshFormat"){
                std::getline(in,line);
                auto v = detail::split_ll(line);
                // line is: <version> <file-type> <data-size>
                // We parse with stringstream to keep decimal version.
                {
                    std::stringstream ss(line);
                    ss >> msh_version;
                    int file_type=0, data_size=0;
                    ss >> file_type >> data_size;
                    is_binary = (file_type==1);
                }
                // consume $EndMeshFormat
                while(std::getline(in,line)){
                    line = detail::trim(line);
                    if(line=="$EndMeshFormat") break;
                }
                break;
            }
            if(line=="$Nodes" || line=="$Elements"){
                // No MeshFormat section; rewind by seeking to beginning and parse normally.
                in.clear();
                in.seekg(0, std::ios::beg);
                break;
            }
        }
        if(is_binary) throw std::runtime_error("MeshP2: binary .msh not supported (export ASCII).");
        if(msh_version >= 4.0){
            throw std::runtime_error(
                "MeshP2: .msh version >= 4 detected. "
                "Export your mesh as 'MSH 2.2 ASCII' in Gmsh (File > Export, options).");
        }

        // reset
        in.clear();
        in.seekg(0, std::ios::beg);

        // Temporary storage with gmsh node tags
        std::unordered_map<long long, int> tag_to_idx; // nodeTag -> index in nodes
        nodes.clear();
        triangles.clear();

        // boundary edges: (v1,v2) -> physical tag (reference)
        std::unordered_map<detail::EdgeKey, int, detail::EdgeKeyHash> boundary_ref;

        // --- Find and parse $Nodes ---
        while(std::getline(in,line)){
            line = detail::trim(line);
            if(line=="$Nodes") break;
        }
        if(!in) throw std::runtime_error("MeshP2: missing $Nodes section");

        std::getline(in,line);
        long long nb_nodes = std::stoll(detail::trim(line));
        nodes.reserve(static_cast<size_t>(nb_nodes));

        for(long long i=0;i<nb_nodes;i++){
            std::getline(in,line);
            auto vals = detail::split_ll(line);
            if(vals.size()<4) throw std::runtime_error("MeshP2: invalid node line: " + line);
            long long tag = vals[0];
            double x=0,y=0;
            {
                std::stringstream ss(line);
                ss >> tag >> x >> y; // ignore z
            }
            Point2D p; p.x=x; p.y=y; p.id=static_cast<int>(nodes.size());
            nodes.push_back(p);
            tag_to_idx[tag]=p.id;
        }
        // consume $EndNodes
        while(std::getline(in,line)){
            line = detail::trim(line);
            if(line=="$EndNodes") break;
        }

        // --- Find and parse $Elements ---
        while(std::getline(in,line)){
            line = detail::trim(line);
            if(line=="$Elements") break;
        }
        if(!in) throw std::runtime_error("MeshP2: missing $Elements section");

        std::getline(in,line);
        long long nb_elems = std::stoll(detail::trim(line));

        struct TriP1Tmp{ int v1,v2,v3; int ref; };
        std::vector<TriP1Tmp> tri_p1;
        tri_p1.reserve(static_cast<size_t>(nb_elems));

        for(long long e=0;e<nb_elems;e++){
            std::getline(in,line);
            if(!in) break;
            auto toks = detail::split_ll(line);
            if(toks.size()<4) continue;
            int elm_type = static_cast<int>(toks[1]);
            int ntags    = static_cast<int>(toks[2]);

            // In MSH2, first tag is physical entity id (reference).
            int physical = 0;
            if(static_cast<int>(toks.size()) >= 3 + ntags + 1){
                if(ntags>0) physical = static_cast<int>(toks[3]);
            }

            // Nodes start at index (3+ntags)
            int start = 3 + ntags;
            if(elm_type==1){ // 2-node line
                if(static_cast<int>(toks.size()) < start+2) continue;
                int a = tag_to_idx.at(toks[start]);
                int b = tag_to_idx.at(toks[start+1]);
                boundary_ref[detail::EdgeKey(a,b)] = physical;
                // mark node refs (optional)
                nodes[a].ref = (nodes[a].ref==0? physical : nodes[a].ref);
                nodes[b].ref = (nodes[b].ref==0? physical : nodes[b].ref);
            } else if(elm_type==2){ // 3-node triangle
                if(static_cast<int>(toks.size()) < start+3) continue;
                int v1 = tag_to_idx.at(toks[start]);
                int v2 = tag_to_idx.at(toks[start+1]);
                int v3 = tag_to_idx.at(toks[start+2]);
                tri_p1.push_back({v1,v2,v3,physical});
            } else {
                // ignore other element types
            }
        }
        // done, no need to parse $EndElements strictly

        // --- Enrich to P2 (create mid-edge nodes shared between triangles) ---
        std::unordered_map<detail::EdgeKey, int, detail::EdgeKeyHash> edge_mid;

        auto get_mid = [&](int a, int b)->int{
            detail::EdgeKey k(a,b);
            auto it = edge_mid.find(k);
            if(it!=edge_mid.end()) return it->second;
            const auto& pa = nodes[a];
            const auto& pb = nodes[b];
            Point2D pm;
            pm.x = 0.5*(pa.x+pb.x);
            pm.y = 0.5*(pa.y+pb.y);
            pm.id = static_cast<int>(nodes.size());
            // If this edge is on boundary, propagate its ref to the midpoint.
            auto itb = boundary_ref.find(k);
            if(itb!=boundary_ref.end()) pm.ref = itb->second;
            nodes.push_back(pm);
            edge_mid[k]=pm.id;
            return pm.id;
        };

        triangles.reserve(tri_p1.size());
        for(const auto& t : tri_p1){
            TriangleP2 T;
            T.ref = t.ref;
            T.node_ids[0]=t.v1;
            T.node_ids[1]=t.v2;
            T.node_ids[2]=t.v3;
            T.node_ids[3]=get_mid(t.v1,t.v2);
            T.node_ids[4]=get_mid(t.v2,t.v3);
            T.node_ids[5]=get_mid(t.v1,t.v3);

            // Boundary refs per edge (vertex edges)
            auto e12 = boundary_ref.find(detail::EdgeKey(t.v1,t.v2));
            auto e23 = boundary_ref.find(detail::EdgeKey(t.v2,t.v3));
            auto e31 = boundary_ref.find(detail::EdgeKey(t.v3,t.v1));
            T.edge_ref[0] = (e12==boundary_ref.end()? 0 : e12->second);
            T.edge_ref[1] = (e23==boundary_ref.end()? 0 : e23->second);
            T.edge_ref[2] = (e31==boundary_ref.end()? 0 : e31->second);

            triangles.push_back(T);
        }

        set_defect_tags(defect_tags);
    }

    // Mark defects from triangle refs
    void set_defect_tags(const std::vector<int>& defect_tags){
        if(defect_tags.empty()){
            for(auto& t: triangles) t.is_defect=false;
            return;
        }
        std::unordered_map<int,bool> is_def;
        for(int r: defect_tags) is_def[r]=true;
        for(auto& t: triangles) t.is_defect = (is_def.find(t.ref)!=is_def.end());
    }

    // Convenience: number of P2 dofs (one per node in this simplified project)
    std::size_t ndof() const { return nodes.size(); }

    // Export to a Matlab/Octave .m script (mesh only).
    // Creates variables:
    //   Coor: (Nbpt x 2)
    //   Tri : (Nbtri x 6)   (1-based indices for Matlab)
    //   RefTri: (Nbtri x 1)
    //   EdgeRef: (Nbtri x 3)
    void write_matlab_mesh_m(const std::string& out_m_file) const {
        std::ofstream out(out_m_file);
        if(!out) throw std::runtime_error("MeshP2: cannot write: " + out_m_file);

        out << "% Auto-generated by MeshP2::write_matlab_mesh_m\n";
        out << "Coor = [ ...\n";
        for(const auto& p: nodes){
            out << p.x << " " << p.y << ";\n";
        }
        out << "];\n\n";

        out << "Tri = [ ...\n";
        for(const auto& t: triangles){
            out << (t.node_ids[0]+1) << " " << (t.node_ids[1]+1) << " " << (t.node_ids[2]+1) << " "
                << (t.node_ids[3]+1) << " " << (t.node_ids[4]+1) << " " << (t.node_ids[5]+1) << ";\n";
        }
        out << "];\n\n";

        out << "RefTri = [ ...\n";
        for(const auto& t: triangles) out << t.ref << ";\n";
        out << "];\n\n";

        out << "EdgeRef = [ ...\n";
        for(const auto& t: triangles){
            out << t.edge_ref[0] << " " << t.edge_ref[1] << " " << t.edge_ref[2] << ";\n";
        }
        out << "];\n\n";

        out << "% Quick visualization (P1 faces only):\n";
        out << "% trisurf(Tri(:,1:3), Coor(:,1), Coor(:,2), 0*Coor(:,1)); view(2); axis equal; shading interp;\n";
        out << "% title('Mesh (P1 view)');\n";
    }

    // Export nodal complex field to Matlab/Octave .m script.
    // Creates variables:
    //   Ure, Uim (Nbpt x 1), and optionally Umag.
    void write_matlab_field_m(const std::string& out_m_file,
                              const std::vector<complexe>& U,
                              const std::string& var_prefix="U") const {
        if(U.size()!=nodes.size())
            throw std::invalid_argument("write_matlab_field_m: U size must equal Nbpt (nodes.size()).");
        std::ofstream out(out_m_file);
        if(!out) throw std::runtime_error("MeshP2: cannot write: " + out_m_file);

        out << "% Auto-generated by MeshP2::write_matlab_field_m\n";
        out << var_prefix << "re = [ ...\n";
        for(const auto& u: U) out << std::real(u) << ";\n";
        out << "];\n\n";
        out << var_prefix << "im = [ ...\n";
        for(const auto& u: U) out << std::imag(u) << ";\n";
        out << "];\n\n";
        out << var_prefix << "mag = sqrt("<<var_prefix<<"re.^2 + "<<var_prefix<<"im.^2);\n";
    }

    // Export a simple ASCII point cloud: x y Re(U) Im(U)
    void write_point_cloud_txt(const std::string& out_txt,
                               const std::vector<complexe>& U) const {
        if(U.size()!=nodes.size())
            throw std::invalid_argument("write_point_cloud_txt: U size must equal Nbpt.");
        std::ofstream out(out_txt);
        if(!out) throw std::runtime_error("MeshP2: cannot write: " + out_txt);
        out << "# x y ReU ImU\n";
        for(std::size_t i=0;i<nodes.size();++i){
            out << nodes[i].x << " " << nodes[i].y << " " << std::real(U[i]) << " " << std::imag(U[i]) << "\n";
        }
    }
};

} // namespace usim

#endif // MESH_HPP
