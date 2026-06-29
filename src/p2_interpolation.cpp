#include "p2_interpolation.hpp"

#include <algorithm>
#include <cmath>

namespace usim {

P2PointLocation locate_point_in_triangle_p2(const MeshP2& mesh, std::size_t triangle_id,
                                            double x, double y, double tolerance) {
    if (triangle_id >= mesh.triangles.size()) return {};
    const auto& triangle = mesh.triangles[triangle_id];
    const auto& p0 = mesh.nodes[triangle.node_ids[0]];
    const auto& p1 = mesh.nodes[triangle.node_ids[1]];
    const auto& p2 = mesh.nodes[triangle.node_ids[2]];

    const double min_x = std::min({p0.x, p1.x, p2.x}) - tolerance;
    const double max_x = std::max({p0.x, p1.x, p2.x}) + tolerance;
    const double min_y = std::min({p0.y, p1.y, p2.y}) - tolerance;
    const double max_y = std::max({p0.y, p1.y, p2.y}) + tolerance;
    if (x < min_x || x > max_x || y < min_y || y > max_y) return {};

    const double j00 = p1.x - p0.x;
    const double j01 = p2.x - p0.x;
    const double j10 = p1.y - p0.y;
    const double j11 = p2.y - p0.y;
    const double determinant = j00 * j11 - j01 * j10;
    if (std::abs(determinant) <= std::numeric_limits<double>::epsilon()) return {};

    const double dx = x - p0.x;
    const double dy = y - p0.y;
    const double xi = (j11 * dx - j01 * dy) / determinant;
    const double eta = (-j10 * dx + j00 * dy) / determinant;
    const double lambda1 = 1.0 - xi - eta;

    if (xi >= -tolerance && eta >= -tolerance && lambda1 >= -tolerance &&
        xi <= 1.0 + tolerance && eta <= 1.0 + tolerance) {
        return {triangle_id, xi, eta};
    }
    return {};
}

P2PointLocation locate_point_p2(const MeshP2& mesh, double x, double y,
                                double tolerance) {
    for (std::size_t triangle_id = 0; triangle_id < mesh.triangles.size(); ++triangle_id) {
        const auto location =
            locate_point_in_triangle_p2(mesh, triangle_id, x, y, tolerance);
        if (location.found()) return location;
    }
    return {};
}

} // namespace usim
