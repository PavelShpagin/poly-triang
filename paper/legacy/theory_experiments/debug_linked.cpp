#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include "variant_fast_linked.hpp"

std::vector<fast_linked::Point> make_random(int n, uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(0.0, 1000.0);
  std::vector<double> angles(n);
  for (int i = 0; i < n; ++i) angles[i] = dist(rng);
  std::sort(angles.begin(), angles.end());
  std::vector<fast_linked::Point> pts(n);
  for (int i = 0; i < n; ++i) {
    double a = angles[i] * 2.0 * M_PI / 1000.0;
    double r = 100.0 + dist(rng) * 0.5;
    pts[i].x = 500.0 + r * std::cos(a);
    pts[i].y = 500.0 + r * std::sin(a);
  }
  return pts;
}

int main() {
  for (int n : {10, 50, 100}) {
    auto pts = make_random(n, 42);
    fast_linked::Triangulator tri;
    tri.triangulate(pts);
    std::cout << "n=" << n << " triangles=" << tri.triangles.size() 
              << " expected=" << n - 2 
              << " reflex=" << tri.reflex_count << std::endl;
  }
  return 0;
}
