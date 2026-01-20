/**
 * Clean O(n + k log k) triangulation CLI.
 * Uses chain_triangulate_clean.hpp - no heuristics, no fallbacks.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <cstring>

#include "chain_triangulate_clean.hpp"

static int count_local_maxima(const std::vector<chain_clean::Point>& pts) {
  const int n = static_cast<int>(pts.size());
  if (n < 3) return 0;
  int k = 0;
  for (int i = 0; i < n; ++i) {
    const int p = (i - 1 + n) % n;
    const int nx = (i + 1) % n;
    auto below = [&](int a, int b) {
      if (pts[a].y < pts[b].y) return true;
      if (pts[a].y > pts[b].y) return false;
      if (pts[a].x > pts[b].x) return true;
      if (pts[a].x < pts[b].x) return false;
      return a > b;
    };
    if (below(p, i) && below(nx, i)) ++k;
  }
  return k;
}

int main(int argc, char* argv[]) {
  std::string input_file, output_file;

  for (int i = 1; i < argc; i++) {
    if ((strcmp(argv[i], "--input") == 0 || strcmp(argv[i], "-i") == 0) && i + 1 < argc) {
      input_file = argv[++i];
    } else if ((strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0) && i + 1 < argc) {
      output_file = argv[++i];
    }
  }

  if (input_file.empty() || output_file.empty()) {
    std::cerr << "Usage: " << argv[0] << " --input <polygon.poly> --output <output.tri>\n";
    return 1;
  }

  std::ifstream fin(input_file);
  if (!fin) {
    std::cerr << "Error: Cannot open " << input_file << "\n";
    return 1;
  }

  int n;
  fin >> n;
  std::vector<chain_clean::Point> pts(n);
  for (int i = 0; i < n; i++) {
    fin >> pts[i].x >> pts[i].y;
    pts[i].index = i;
  }
  fin.close();

  int k_count = count_local_maxima(pts);

  chain_clean::Triangulator tri;
  auto start = std::chrono::high_resolution_clock::now();
  try {
    tri.triangulate(pts);
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 2;
  }
  auto end = std::chrono::high_resolution_clock::now();

  double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();

  std::ofstream fout(output_file);
  if (!fout) {
    std::cerr << "Error: Cannot write to " << output_file << "\n";
    return 1;
  }
  fout << "# vertices\n" << n << "\n";
  for (int i = 0; i < n; i++) {
    fout << pts[i].x << " " << pts[i].y << "\n";
  }
  fout << "# triangles\n" << tri.triangles.size() << "\n";
  for (const auto& t : tri.triangles) {
    fout << t.v0 << " " << t.v1 << " " << t.v2 << "\n";
  }
  fout.close();

  std::cout << "chain_clean,vertices=" << n
            << ",triangles=" << tri.triangles.size()
            << ",expected=" << (n - 2)
            << ",k_count=" << k_count
            << ",reflex_count=" << tri.reflex_count
            << ",diagonals=" << tri.diagonal_count
            << ",faces=" << tri.debug_face_count
            << ",interior_faces=" << tri.debug_interior_faces
            << ",time_ms=" << elapsed_ms << "\n";

  return (static_cast<int>(tri.triangles.size()) == n - 2) ? 0 : 2;
}
