#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include "../../methods/include/reflex_fast.hpp"

std::vector<reflex_fast::Point> make_random(int n, uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(0.0, 1000.0);
  std::vector<double> angles(n);
  for (int i = 0; i < n; ++i) angles[i] = dist(rng);
  std::sort(angles.begin(), angles.end());
  std::vector<reflex_fast::Point> pts(n);
  for (int i = 0; i < n; ++i) {
    double a = angles[i] * 2.0 * M_PI / 1000.0;
    double r = 100.0 + dist(rng) * 0.5;
    pts[i].x = 500.0 + r * std::cos(a);
    pts[i].y = 500.0 + r * std::sin(a);
    pts[i].index = i;
  }
  return pts;
}

int main() {
  for (int n : {10, 20, 50, 100, 500, 1000}) {
    int ok = 0, fail = 0;
    for (int seed = 0; seed < 20; ++seed) {
      auto pts = make_random(n, seed);
      reflex_fast::FastTriangulator tri;
      tri.triangulate(pts);
      if (static_cast<int>(tri.triangles.size()) == n - 2) ok++;
      else fail++;
    }
    std::cout << "n=" << n << " ok=" << ok << "/20 fail=" << fail << std::endl;
  }
  return 0;
}

