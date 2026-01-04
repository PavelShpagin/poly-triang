#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include "variant_fast_linked.hpp"
#include "../baselines/polypartition/src/polypartition.h"

using Clock = std::chrono::high_resolution_clock;

static std::vector<fast_linked::Point> make_random(int n, uint32_t seed) {
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
  constexpr int N = 10000;
  constexpr int TRIALS = 100;
  
  double sum_ours = 0, sum_pp = 0;
  int ours_correct = 0, pp_correct = 0, both_correct = 0;
  int ours_wins = 0;
  
  for (int t = 0; t < TRIALS; ++t) {
    auto pts = make_random(N, 1000 + t);
    
    // PolyPartition first (to check if it succeeds)
    TPPLPoly poly;
    poly.Init(N);
    for (int i = 0; i < N; ++i) {
      poly[i].x = pts[i].x;
      poly[i].y = pts[i].y;
    }
    poly.SetOrientation(TPPL_ORIENTATION_CCW);
    std::list<TPPLPoly> polys, triangles_pp;
    polys.push_back(poly);
    TPPLPartition pp;
    auto t0 = Clock::now();
    int pp_result = pp.Triangulate_MONO(&polys, &triangles_pp);
    auto t1 = Clock::now();
    double pp_time = std::chrono::duration<double, std::milli>(t1 - t0).count();
    bool pp_ok = (pp_result == 1 && triangles_pp.size() == N - 2);
    if (pp_ok) pp_correct++;
    
    // Our algorithm
    fast_linked::Triangulator tri;
    t0 = Clock::now();
    tri.triangulate(pts);
    t1 = Clock::now();
    double ours_time = std::chrono::duration<double, std::milli>(t1 - t0).count();
    bool ours_ok = (tri.triangles.size() == N - 2);
    if (ours_ok) ours_correct++;
    
    if (pp_ok && ours_ok) {
      both_correct++;
      sum_pp += pp_time;
      sum_ours += ours_time;
      if (ours_time < pp_time) ours_wins++;
    }
  }
  
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Fair Benchmark: n=" << N << ", " << TRIALS << " trials\n";
  std::cout << "============================================\n\n";
  
  std::cout << "Correctness:\n";
  std::cout << "  PolyPartition:  " << pp_correct << "/" << TRIALS << " (" 
            << (100.0 * pp_correct / TRIALS) << "%)\n";
  std::cout << "  Ours:           " << ours_correct << "/" << TRIALS << " ("
            << (100.0 * ours_correct / TRIALS) << "%)\n";
  std::cout << "  Both correct:   " << both_correct << "/" << TRIALS << "\n\n";
  
  if (both_correct > 0) {
    std::cout << "Performance (on " << both_correct << " both-correct cases):\n";
    std::cout << "  PolyPartition:  " << (sum_pp / both_correct) << " ms avg\n";
    std::cout << "  Ours:           " << (sum_ours / both_correct) << " ms avg\n";
    std::cout << "  Ratio:          " << (sum_ours / sum_pp) << "x\n";
    std::cout << "  Ours wins:      " << ours_wins << "/" << both_correct << " ("
              << (100.0 * ours_wins / both_correct) << "%)\n";
    
    if (sum_ours < sum_pp) {
      double speedup = 100.0 * (1.0 - sum_ours / sum_pp);
      std::cout << "\n  *** WE ARE " << speedup << "% FASTER! ***\n";
    }
  }
  
  return 0;
}

