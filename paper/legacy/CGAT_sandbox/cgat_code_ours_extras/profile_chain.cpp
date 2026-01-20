/**
 * Profiling CLI: prints phase breakdown for chain-only triangulation.
 *
 * Build:
 *   g++ -O3 -std=c++17 -DREFLEX_TRI_TIMING profile_chain.cpp -o profile_chain
 */

#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>

#include "reflex_chain_triangulate.hpp"

int main(int argc, char* argv[]) {
  std::string input;
  for (int i = 1; i < argc; ++i) {
    if ((std::strcmp(argv[i], "--input") == 0 || std::strcmp(argv[i], "-i") == 0) && i + 1 < argc) {
      input = argv[++i];
    }
  }
  if (input.empty()) {
    std::cerr << "Usage: " << argv[0] << " --input <poly.poly>\n";
    return 1;
  }

  std::ifstream fin(input);
  if (!fin) {
    std::cerr << "Error: cannot open " << input << "\n";
    return 1;
  }

  int n = 0;
  fin >> n;
  std::vector<reflex_tri::Point> pts(n);
  for (int i = 0; i < n; ++i) {
    fin >> pts[i].x >> pts[i].y;
    pts[i].index = i;
  }

  reflex_tri::Triangulator tri;
  auto t0 = std::chrono::high_resolution_clock::now();
  try {
    tri.triangulate(pts);
  } catch (const std::exception& e) {
    std::cerr << "triangulate() threw: " << e.what() << "\n";
  }
  auto t1 = std::chrono::high_resolution_clock::now();

  const double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

  std::cout << "n=" << n
            << " triangles=" << tri.debug_triangles().size()
            << " total_ms=" << total_ms
            << " phase1_ms=" << tri.phase1_ms()
            << " phase2_ms=" << tri.phase2_ms()
            << " phase3_ms=" << tri.phase3_ms()
            << "\n";

  return 0;
}

