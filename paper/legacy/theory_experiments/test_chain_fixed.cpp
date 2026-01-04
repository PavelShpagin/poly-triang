#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <chrono>
#include "reflex_chain_v2.hpp"
#include "../baselines/polypartition/src/polypartition.h"
#include "../baselines/polypartition/src/polypartition.cpp"

std::vector<reflex_chain_v2::Point> make_random(int n, uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(0.0, 1000.0);
  std::vector<double> angles(n);
  for (int i = 0; i < n; ++i) angles[i] = dist(rng);
  std::sort(angles.begin(), angles.end());
  std::vector<reflex_chain_v2::Point> pts(n);
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
  std::cout << "Testing O(n + r log r) chain-based implementation v2\n\n";
  
  std::cout << "=== Correctness ===\n";
  for (int n : {10, 20, 50, 100, 500, 1000, 5000, 10000}) {
    int ok = 0, trials = (n <= 100) ? 50 : 20;
    for (int seed = 0; seed < trials; ++seed) {
      auto pts = make_random(n, seed);
      reflex_chain_v2::Triangulator tri;
      tri.triangulate(pts);
      if (static_cast<int>(tri.triangles.size()) == n - 2) ok++;
    }
    std::cout << "n=" << n << " ok=" << ok << "/" << trials << "\n";
  }
  
  // Performance comparison
  std::cout << "\n=== Performance (n=10000, 50 trials) ===\n";
  int n = 10000;
  double our_total = 0, pp_total = 0;
  int our_ok = 0, pp_ok = 0, both_ok = 0;
  int our_wins = 0;
  
  for (int seed = 0; seed < 50; ++seed) {
    auto pts = make_random(n, seed);
    
    // Our algorithm
    reflex_chain_v2::Triangulator tri;
    auto t1 = std::chrono::high_resolution_clock::now();
    tri.triangulate(pts);
    auto t2 = std::chrono::high_resolution_clock::now();
    double our_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    bool our_correct = (static_cast<int>(tri.triangles.size()) == n - 2);
    
    // PolyPartition
    TPPLPoly poly;
    poly.Init(n);
    for (int i = 0; i < n; ++i) {
      poly[i].x = pts[i].x;
      poly[i].y = pts[i].y;
    }
    poly.SetOrientation(TPPL_ORIENTATION_CCW);
    
    std::list<TPPLPoly> result;
    TPPLPartition part;
    t1 = std::chrono::high_resolution_clock::now();
    int ret = part.Triangulate_MONO(&poly, &result);
    t2 = std::chrono::high_resolution_clock::now();
    double pp_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    bool pp_correct = ret && (static_cast<int>(result.size()) == n - 2);
    
    if (our_correct) our_ok++;
    if (pp_correct) pp_ok++;
    
    if (our_correct && pp_correct) {
      both_ok++;
      our_total += our_ms;
      pp_total += pp_ms;
      if (our_ms < pp_ms) our_wins++;
    }
  }
  
  std::cout << "Ours:          " << our_ok << "/50 correct\n";
  std::cout << "PolyPartition: " << pp_ok << "/50 correct\n";
  std::cout << "Both correct:  " << both_ok << "/50\n";
  
  if (both_ok > 0) {
    std::cout << "\nPerformance (on " << both_ok << " both-correct cases):\n";
    std::cout << "  Ours:          " << (our_total / both_ok) << " ms avg\n";
    std::cout << "  PolyPartition: " << (pp_total / both_ok) << " ms avg\n";
    std::cout << "  Ratio:         " << (our_total / pp_total) << "x\n";
    std::cout << "  Ours wins:     " << our_wins << "/" << both_ok << "\n";
    double speedup = 100.0 * (1.0 - (our_total / both_ok) / (pp_total / both_ok));
    if (speedup > 0) {
      std::cout << "\n  *** WE ARE " << speedup << "% FASTER! ***\n";
    } else {
      std::cout << "\n  *** We are " << -speedup << "% slower ***\n";
    }
  }
  
  return 0;
}

