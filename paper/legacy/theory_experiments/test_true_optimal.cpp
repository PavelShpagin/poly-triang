#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <chrono>
#include "variant_true_optimal.hpp"
#include "../baselines/polypartition/src/polypartition.h"
#include "../baselines/polypartition/src/polypartition.cpp"

std::vector<true_optimal::Point> make_random(int n, uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(0.0, 1000.0);
  std::vector<double> angles(n);
  for (int i = 0; i < n; ++i) angles[i] = dist(rng);
  std::sort(angles.begin(), angles.end());
  std::vector<true_optimal::Point> pts(n);
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
  std::cout << "Testing TRUE O(n + r log r) implementation\n\n";
  
  // Correctness test
  std::cout << "=== Correctness ===\n";
  for (int n : {10, 20, 50, 100, 500, 1000}) {
    int ok = 0;
    for (int seed = 0; seed < 20; ++seed) {
      auto pts = make_random(n, seed);
      true_optimal::Triangulator tri;
      tri.triangulate(pts);
      if (static_cast<int>(tri.triangles.size()) == n - 2) ok++;
    }
    std::cout << "n=" << n << " ok=" << ok << "/20\n";
  }
  
  // Performance comparison
  std::cout << "\n=== Performance (n=10000, 50 trials) ===\n";
  int n = 10000;
  double our_total = 0, pp_total = 0;
  int our_ok = 0, pp_ok = 0;
  
  for (int seed = 0; seed < 50; ++seed) {
    auto pts = make_random(n, seed);
    
    // Our algorithm
    true_optimal::Triangulator tri;
    auto t1 = std::chrono::high_resolution_clock::now();
    tri.triangulate(pts);
    auto t2 = std::chrono::high_resolution_clock::now();
    double our_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    
    if (static_cast<int>(tri.triangles.size()) == n - 2) {
      our_total += our_ms;
      our_ok++;
    }
    
    // PolyPartition
    TPPLPoly poly;
    poly.Init(n);
    for (int i = 0; i < n; ++i) {
      poly[i].x = pts[i].x;
      poly[i].y = pts[i].y;
    }
    poly.SetOrientation(TPPL_CCW);
    
    std::list<TPPLPoly> result;
    TPPLPartition part;
    t1 = std::chrono::high_resolution_clock::now();
    int ret = part.Triangulate_MONO(&poly, &result);
    t2 = std::chrono::high_resolution_clock::now();
    double pp_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    
    if (ret && static_cast<int>(result.size()) == n - 2) {
      pp_total += pp_ms;
      pp_ok++;
    }
  }
  
  std::cout << "Ours:         " << our_ok << "/50 correct, avg " << (our_ok > 0 ? our_total/our_ok : 0) << " ms\n";
  std::cout << "PolyPartition: " << pp_ok << "/50 correct, avg " << (pp_ok > 0 ? pp_total/pp_ok : 0) << " ms\n";
  if (our_ok > 0 && pp_ok > 0) {
    double speedup = 100.0 * (1.0 - (our_total/our_ok) / (pp_total/pp_ok));
    std::cout << "Speedup: " << speedup << "%\n";
  }
  
  return 0;
}

