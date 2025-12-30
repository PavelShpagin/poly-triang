#define POLYTRI_IMPLEMENTATION
#include <algorithm>
#include <iostream>
#include <vector>

#include "polygon_io.hpp"
#include "polytri/polytri.hpp"

static double signed_area(const baselines::Polygon& polygon) {
  double area = 0.0;
  const std::size_t n = polygon.size();
  for (std::size_t i = 0; i < n; ++i) {
    const auto& a = polygon[i];
    const auto& b = polygon[(i + 1) % n];
    area += a.x * b.y - b.x * a.y;
  }
  return area * 0.5;
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: polytri_read <file.poly>\n";
    return 2;
  }
  auto polygon = baselines::read_polygon(argv[1]);
  if (signed_area(polygon) < 0.0) {
    std::reverse(polygon.begin(), polygon.end());
  }

  std::vector<std::vector<PolyTri::vertex_t>> contours(1);
  contours[0].reserve(polygon.size());
  for (const auto& p : polygon) {
    contours[0].emplace_back(p.x, p.y);
  }

  auto indices = PolyTri::Triangulate(contours);
  std::cout << "ok indices=" << indices.size() << "\n";
  return 0;
}


