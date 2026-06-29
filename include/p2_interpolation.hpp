#ifndef P2_INTERPOLATION_HPP
#define P2_INTERPOLATION_HPP

#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "mesh.hpp"

namespace usim {

struct P2PointLocation {
    std::size_t triangle_id = std::numeric_limits<std::size_t>::max();
    double xi = 0.0;
    double eta = 0.0;

    bool found() const {
        return triangle_id != std::numeric_limits<std::size_t>::max();
    }
};

// Locate a physical point in the affine P1 geometry underlying a P2 triangle.
P2PointLocation locate_point_in_triangle_p2(const MeshP2& mesh, std::size_t triangle_id,
                                            double x, double y,
                                            double tolerance = 1e-10);

P2PointLocation locate_point_p2(const MeshP2& mesh, double x, double y,
                                double tolerance = 1e-10);

inline std::array<double, 6> p2_shape_functions(double xi, double eta) {
    const double lambda1 = 1.0 - xi - eta;
    const double lambda2 = xi;
    const double lambda3 = eta;
    return {
        lambda1 * (2.0 * lambda1 - 1.0),
        lambda2 * (2.0 * lambda2 - 1.0),
        lambda3 * (2.0 * lambda3 - 1.0),
        4.0 * lambda1 * lambda2,
        4.0 * lambda2 * lambda3,
        4.0 * lambda1 * lambda3
    };
}

template <typename T>
T interpolate_p2(const MeshP2& mesh, const std::vector<T>& nodal_values,
                 const P2PointLocation& location) {
    if (!location.found() || location.triangle_id >= mesh.triangles.size()) {
        throw std::invalid_argument("interpolate_p2: invalid point location");
    }
    if (nodal_values.size() != mesh.ndof()) {
        throw std::invalid_argument("interpolate_p2: nodal_values size must equal mesh.ndof()");
    }

    const auto phi = p2_shape_functions(location.xi, location.eta);
    const auto& triangle = mesh.triangles[location.triangle_id];
    T value{};
    for (std::size_t i = 0; i < phi.size(); ++i) {
        value += nodal_values[triangle.node_ids[i]] * phi[i];
    }
    return value;
}

} // namespace usim

#endif // P2_INTERPOLATION_HPP
