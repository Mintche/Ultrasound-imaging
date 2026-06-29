#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "p2_interpolation.hpp"

using namespace usim;

int main() {
    MeshP2 mesh;
    mesh.nodes = {
        {0.0, 0.0, 0, 0}, {1.0, 0.0, 1, 0}, {0.0, 1.0, 2, 0},
        {0.5, 0.0, 3, 0}, {0.5, 0.5, 4, 0}, {0.0, 0.5, 5, 0}
    };
    TriangleP2 triangle;
    triangle.node_ids = {0, 1, 2, 3, 4, 5};
    mesh.triangles.push_back(triangle);

    // Polynomial P2 reproduced exactly by the six nodal values.
    auto polynomial = [](double x, double y) {
        return std::complex<double>(1.0 + 2.0*x - 3.0*y + 4.0*x*x + 2.0*x*y,
                                    -2.0 + x + y*y);
    };
    std::vector<complexe> values;
    for (const auto& node : mesh.nodes) values.push_back(polynomial(node.x, node.y));

    const double x = 0.2;
    const double y = 0.3;
    const auto location = locate_point_p2(mesh, x, y);
    assert(location.found());
    assert(std::abs(interpolate_p2(mesh, values, location) - polynomial(x, y)) < 1e-12);
    assert(!locate_point_p2(mesh, 0.8, 0.8).found());

    std::cout << "[OK] Interpolation P2 exacte pour un polynome quadratique.\n";
    return 0;
}
