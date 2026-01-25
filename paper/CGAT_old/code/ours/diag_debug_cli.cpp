/**
 * Debug tool: dump diagonal provenance for the chain-based sweep.
 *
 * Usage:
 *   ./diag_debug_cli --input poly.poly
 *   ./diag_debug_cli --input poly.poly --filter 15,18
 *
 * Prints diag records and pending records to stdout.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>

#include "reflex_chain_triangulate.hpp"

static const char* vtype_name(int t) {
  // Matches enum order in reflex_chain_triangulate.hpp: Start, End, Split, Merge, Regular
  switch (t) {
    case 0: return "Start";
    case 1: return "End";
    case 2: return "Split";
    case 3: return "Merge";
    case 4: return "Regular";
    default: return "Unknown";
  }
}

static bool parse_filter(const std::string& s, int& a, int& b) {
  const auto pos = s.find(',');
  if (pos == std::string::npos) return false;
  try {
    a = std::stoi(s.substr(0, pos));
    b = std::stoi(s.substr(pos + 1));
  } catch (...) {
    return false;
  }
  return true;
}

int main(int argc, char* argv[]) {
  std::string input_file;
  bool do_filter = false;
  int fa = -1, fb = -1;

  for (int i = 1; i < argc; ++i) {
    if ((std::strcmp(argv[i], "--input") == 0 || std::strcmp(argv[i], "-i") == 0) && i + 1 < argc) {
      input_file = argv[++i];
    } else if ((std::strcmp(argv[i], "--filter") == 0) && i + 1 < argc) {
      do_filter = parse_filter(argv[++i], fa, fb);
    }
  }

  if (input_file.empty()) {
    std::cerr << "Usage: " << argv[0] << " --input <poly.poly> [--filter a,b]\n";
    return 1;
  }

  std::ifstream fin(input_file);
  if (!fin) {
    std::cerr << "Error: cannot open " << input_file << "\n";
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
  try {
    (void)tri.triangulate(pts);
  } catch (const std::exception& e) {
    std::cerr << "triangulate() threw: " << e.what() << "\n";
  }

  std::cout << "n=" << n
            << " reflex=" << tri.reflex_count()
            << " diagonals=" << tri.debug_diagonals().size()
            << " diag_records=" << tri.debug_diag_records().size()
            << " pending_records=" << tri.debug_pending_records().size()
            << "\n";

  std::cout << "\n# diag_records: a,b,event_v,event_type,chosen_chain,target_v,reason\n";
  for (const auto& r : tri.debug_diag_records()) {
    if (do_filter) {
      const int a = r.a, b = r.b;
      if (!((a == fa && b == fb) || (a == fb && b == fa))) continue;
    }
    std::cout << r.a << "," << r.b
              << "," << r.event_v
              << "," << vtype_name(r.event_type)
              << "," << r.chosen_chain
              << "," << r.target_v
              << "," << r.reason
              << "\n";
  }

  // Always print the raw diagonal set (independent of NDEBUG), so correctness
  // checkers can validate non-crossing without needing debug provenance.
  std::cout << "\n# diagonals: a,b\n";
  for (const auto& d : tri.debug_diagonals()) {
    std::cout << d.first << "," << d.second << "\n";
  }

  std::cout << "\n# pending_records: chain_id,pending_v,event_v,reason\n";
  for (const auto& r : tri.debug_pending_records()) {
    std::cout << r.chain_id
              << "," << r.pending_v
              << "," << r.event_v
              << "," << r.reason
              << "\n";
  }

  return 0;
}

