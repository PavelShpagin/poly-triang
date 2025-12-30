#define POLYTRI_IMPLEMENTATION
#include <iostream>
#include <vector>

#include "polytri/polytri.hpp"

int main() {
  std::vector<std::vector<PolyTri::vertex_t>> polys = {
      {{-10.0, -9.0}, {11.0, -12.0}, {0.0, 8.0}, {-5.0, 11.0}}};
  auto indices = PolyTri::Triangulate(polys);
  std::cout << "indices=" << indices.size() << "\n";
  for (auto i : indices) {
    std::cout << i << " ";
  }
  std::cout << "\n";
  return 0;
}


