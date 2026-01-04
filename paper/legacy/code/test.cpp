#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "reflex_triangulate.hpp"

static std::vector<reflex::Point> read_poly(const std::string& path) {
  std::ifstream fin(path);
  if (!fin) throw std::runtime_error("cannot open input: " + path);
  int n = 0;
  fin >> n;
  if (n < 3) throw std::runtime_error("bad n");
  std::vector<reflex::Point> pts(n);
  for (int i = 0; i < n; ++i) {
    fin >> pts[i].x >> pts[i].y;
    pts[i].index = i;
  }
  return pts;
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <polygon.poly>\n";
    return 2;
  }

  auto pts = read_poly(argv[1]);
  const int n = static_cast<int>(pts.size());

  reflex::Triangulator tri;
  tri.triangulate(pts);

  std::cout << "n=" << n
            << " triangles=" << tri.triangles.size()
            << " expected=" << (n - 2)
            << " reflex=" << tri.reflex_count
            << "\n";

  if (static_cast<int>(tri.triangles.size()) != n - 2) return 1;
  return 0;
}
