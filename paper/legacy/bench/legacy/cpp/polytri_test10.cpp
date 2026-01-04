#define POLYTRI_IMPLEMENTATION
#include <cmath>
#include <iostream>
#include <vector>

#include "polytri/polytri.hpp"

int main() {
  const int n = 10;
  std::vector<PolyTri::vertex_t> contour;
  contour.reserve(n);
  for (int i = 0; i < n; ++i) {
    const double a = 2.0 * M_PI * i / n;
    const double r = 100.0 + 0.1 * i;
    contour.emplace_back(r * std::cos(a), r * std::sin(a));
  }

  std::vector<std::vector<PolyTri::vertex_t>> polys(1);
  polys[0] = contour;

  auto indices = PolyTri::Triangulate(polys);
  std::cout << "ok indices=" << indices.size() << "\n";
  return 0;
}


