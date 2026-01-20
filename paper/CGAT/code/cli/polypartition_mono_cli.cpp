// PolyPartition MONO triangulation baseline CLI.
//
// Implements O(n log n) triangulation via monotone partitioning:
//   TPPLPartition::Triangulate_MONO
//
// Input:  .poly (n then n lines x y), simple polygon without holes.
// Output: .tri and one stdout line:
//   garey,vertices=n,triangles=t,time_ms=...
//
// We label it "garey" because Triangulate_MONO corresponds to the classical
// monotone decomposition + triangulation pipeline (Garey et al.).

#include <array>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <time.h>

#include "polypartition.h"

static double cpu_now_ms() {
  struct timespec ts;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
  return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
}

namespace {

struct Args {
  std::string input;
  std::string output;
};

Args parse_args(int argc, char** argv) {
  Args args;
  for (int i = 1; i < argc; ++i) {
    std::string t = argv[i];
    if ((t == "-i" || t == "--input") && i + 1 < argc) {
      args.input = argv[++i];
    } else if ((t == "-o" || t == "--output") && i + 1 < argc) {
      args.output = argv[++i];
    } else {
      throw std::runtime_error("Unknown or incomplete argument: " + t);
    }
  }
  if (args.input.empty()) throw std::runtime_error("Missing --input");
  if (args.output.empty()) throw std::runtime_error("Missing --output");
  return args;
}

struct Pt {
  double x;
  double y;
};

std::vector<Pt> read_poly(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Failed to open input file: " + path);
  int n = 0;
  in >> n;
  if (n < 3) throw std::runtime_error("Polygon must have at least 3 vertices");
  std::vector<Pt> pts;
  pts.reserve(n);
  for (int i = 0; i < n; ++i) {
    Pt p{};
    if (!(in >> p.x >> p.y)) throw std::runtime_error("Malformed polygon file");
    pts.push_back(p);
  }
  return pts;
}

void write_tri(const std::string& path,
               const std::vector<Pt>& poly,
               const std::vector<std::array<int, 3>>& tris) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Failed to open output file: " + path);
  out << "# vertices\n";
  out << poly.size() << "\n";
  out.setf(std::ios::fixed);
  out.precision(10);
  for (const auto& p : poly) {
    out << p.x << " " << p.y << "\n";
  }
  out << "# triangles\n";
  out << tris.size() << "\n";
  for (const auto& t : tris) {
    out << t[0] << " " << t[1] << " " << t[2] << "\n";
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto args = parse_args(argc, argv);
    const auto pts = read_poly(args.input);

    // Time the full computation needed to obtain triangle indices (excluding output I/O),
    // including conversion of the input polygon into the library's representation.
    // This matches the semantics of our "ours" timing (the full CPU time of triangulate()).
    const double t0 = cpu_now_ms();

    TPPLPoly poly;
    poly.Init(static_cast<long>(pts.size()));
    for (long i = 0; i < poly.GetNumPoints(); ++i) {
      poly[i].x = pts[static_cast<std::size_t>(i)].x;
      poly[i].y = pts[static_cast<std::size_t>(i)].y;
      poly[i].id = static_cast<int>(i);  // preserve original index
    }
    poly.SetHole(false);
    poly.SetOrientation(TPPL_ORIENTATION_CCW);

    TPPLPartition partition;
    TPPLPolyList triangles;

    const int ok = partition.Triangulate_MONO(&poly, &triangles);

    if (ok == 0) {
      throw std::runtime_error("Triangulate_MONO failed");
    }

    std::vector<std::array<int, 3>> out_tris;
    out_tris.reserve(triangles.size());
    for (const auto& tri : triangles) {
      if (tri.GetNumPoints() != 3) continue;
      out_tris.push_back({tri[0].id, tri[1].id, tri[2].id});
    }

    const double t1 = cpu_now_ms();
    const double elapsed_ms = (t1 - t0);

    write_tri(args.output, pts, out_tris);
    std::cout << "garey,vertices=" << pts.size()
              << ",triangles=" << out_tris.size()
              << ",time_ms=" << elapsed_ms << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }
}

