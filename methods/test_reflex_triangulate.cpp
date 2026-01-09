// Quick correctness sanity-check for the textbook-style sweep implementation in
// `methods/include/reflex_triangulate.hpp`.
//
// Build & run (WSL):
//   g++ -O3 -std=c++17 -Imethods/include methods/test_reflex_triangulate.cpp -o /tmp/test_reflex_tri
//   /tmp/test_reflex_tri polygons/benchmark/random_500.poly
//
#include <fstream>
#include <iostream>
#include <vector>

#include "reflex_triangulate.hpp"

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <polygon.poly>\n";
    return 2;
  }

  std::ifstream fin(argv[1]);
  if (!fin) {
    std::cerr << "Cannot open: " << argv[1] << "\n";
    return 2;
  }

  int n = 0;
  fin >> n;
  if (n < 3) {
    std::cerr << "Bad n\n";
    return 2;
  }

  std::vector<reflex_tri::Point> pts(n);
  for (int i = 0; i < n; ++i) {
    fin >> pts[i].x >> pts[i].y;
    pts[i].index = i;
  }

  reflex_tri::Triangulator tri;
  const auto& tris = tri.triangulate(pts);

  std::cout << "n=" << n
            << " triangles=" << tris.size()
            << " expected=" << (n - 2)
            << " reflex=" << tri.reflex_count()
            << "\n";
  return (static_cast<int>(tris.size()) == n - 2) ? 0 : 1;
}


