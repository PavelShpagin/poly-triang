#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include "variant_fast_linked.hpp"
#include "../baselines/polypartition/src/polypartition.h"

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
  for (int n : {10, 20, 30, 40, 50, 60, 70, 80, 90, 100}) {
    auto pts = make_random(n, 42);
    
    TPPLPoly poly;
    poly.Init(n);
    for (int i = 0; i < n; ++i) {
      poly[i].x = pts[i].x;
      poly[i].y = pts[i].y;
    }
    poly.SetOrientation(TPPL_ORIENTATION_CCW);
    
    std::list<TPPLPoly> polys, monotone;
    polys.push_back(poly);
    TPPLPartition pp;
    pp.MonotonePartition(&polys, &monotone);
    
    fast_linked::Triangulator tri;
    tri.triangulate(pts);
    
    int expected = n - 2;
    int actual = static_cast<int>(tri.triangles.size());
    bool ok = (expected == actual);
    
    std::cout << "n=" << n 
              << " diags=" << tri.num_diagonals << "/" << monotone.size()-1
              << " faces=" << tri.num_faces << "/" << monotone.size()
              << " tris=" << actual << "/" << expected
              << (ok ? " OK" : " FAIL")
              << std::endl;
  }
  
  return 0;
}

