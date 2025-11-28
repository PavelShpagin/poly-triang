#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace baselines {

struct Point {
    double x;
    double y;
};

using Polygon = std::vector<Point>;

struct Triangle {
    uint32_t a;
    uint32_t b;
    uint32_t c;
};

using Triangles = std::vector<Triangle>;

struct TriangulationResult {
    Triangles triangles;
    double elapsed_ms;
};

Polygon read_polygon(const std::string& path);
void write_triangulation(const Polygon& poly,
                         const Triangles& tris,
                         const std::string& path);

} // namespace baselines

