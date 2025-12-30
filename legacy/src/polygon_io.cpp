#include "polygon_io.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace baselines {

Polygon read_polygon(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open polygon file: " + path);
    }

    size_t n = 0;
    in >> n;
    if (n < 3) {
        throw std::runtime_error("Polygon must have at least 3 vertices");
    }

    Polygon poly;
    poly.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        double x, y;
        if (!(in >> x >> y)) {
            throw std::runtime_error("Malformed polygon file: " + path);
        }
        poly.push_back({x, y});
    }

    return poly;
}

void write_triangulation(const Polygon& poly,
                         const Triangles& tris,
                         const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + path);
    }

    out << "# vertices\n";
    out << poly.size() << "\n";
    out << std::setprecision(10);
    for (const auto& p : poly) {
        out << p.x << " " << p.y << "\n";
    }

    out << "# triangles\n";
    out << tris.size() << "\n";
    for (const auto& t : tris) {
        out << t.a << " " << t.b << " " << t.c << "\n";
    }
}

} // namespace baselines

