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
  int n = 100;
  auto pts = make_random(n, 42);
  
  // Check PolyPartition monotone partition
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
  
  std::cout << "PolyPartition monotone pieces: " << monotone.size() << std::endl;
  int pp_total_verts = 0;
  for (auto& m : monotone) pp_total_verts += m.GetNumPoints();
  std::cout << "PolyPartition total verts in pieces: " << pp_total_verts << std::endl;
  
  // Check with a smaller polygon
  n = 20;
  pts = make_random(n, 42);
  
  poly.Init(n);
  for (int i = 0; i < n; ++i) {
    poly[i].x = pts[i].x;
    poly[i].y = pts[i].y;
  }
  poly.SetOrientation(TPPL_ORIENTATION_CCW);
  
  polys.clear();
  monotone.clear();
  polys.push_back(poly);
  pp.MonotonePartition(&polys, &monotone);
  
  std::cout << "\nSmaller test (n=" << n << "):" << std::endl;
  std::cout << "PolyPartition monotone pieces: " << monotone.size() << std::endl;
  
  // Our algorithm
  fast_linked::Triangulator tri;
  tri.triangulate(pts);
  std::cout << "Our diagonals: " << tri.num_diagonals << std::endl;
  std::cout << "Our faces: " << tri.num_faces << " (expected " << monotone.size() << ")" << std::endl;
  std::cout << "Our triangles: " << tri.triangles.size() << " (expected " << n-2 << ")" << std::endl;
  
  return 0;
}

