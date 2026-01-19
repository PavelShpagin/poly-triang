// CGAL Hertel--Mehlhorn (approximate convex partition) baseline CLI.
//
// Uses CGAL::approx_convex_partition_2 (Partition_2 package), then triangulates
// each convex piece with a fan (O(total vertices)).
//
// Input:  .poly (n then n lines x y), simple polygon without holes.
// Output: .tri and one stdout line:
//   hertel,vertices=n,triangles=t,time_ms=...
//
// NOTE: CGAL's Partition_2 package is GPL-3.0-or-later OR commercial licensed.
// Ensure this baseline is acceptable for your distribution / artifact policy.

#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <time.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/partition_2.h>

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
    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point_2 = K::Point_2;
    // CGAL's approx_convex_partition_2 returns sub-polygons stored in a
    // list-based Polygon_2 type (see Partition_2 internals). Using std::list as
    // the polygon container here avoids type-mismatch issues on output.
    using Polygon_2 = CGAL::Polygon_2<K, std::list<Point_2>>;

    const auto args = parse_args(argc, argv);
    const auto pts = read_poly(args.input);
    const int n = static_cast<int>(pts.size());

    // Coordinate -> original vertex index mapping.
    // We rely on CGAL returning only original vertices in partition polygons.
    std::map<std::pair<double, double>, int> idx_by_xy;
    idx_by_xy.clear();
    for (int i = 0; i < n; ++i) {
      idx_by_xy[{pts[static_cast<std::size_t>(i)].x, pts[static_cast<std::size_t>(i)].y}] = i;
    }

    Polygon_2 poly;
    for (const auto& p : pts) {
      poly.push_back(Point_2(p.x, p.y));
    }
    if (poly.orientation() == CGAL::CLOCKWISE) {
      poly.reverse_orientation();
    }

    std::vector<std::array<int, 3>> out_tris;
    out_tris.reserve(static_cast<std::size_t>(n > 2 ? (n - 2) : 0));

    const double t0 = cpu_now_ms();

    std::list<Polygon_2> parts;
    CGAL::approx_convex_partition_2(poly.vertices_begin(), poly.vertices_end(),
                                   std::back_inserter(parts));

    // Triangulate each convex part by a fan around vertex 0.
    for (const auto& part : parts) {
      const int m = static_cast<int>(part.size());
      if (m < 3) continue;

      std::vector<int> vids;
      vids.reserve(static_cast<std::size_t>(m));
      for (auto vit = part.vertices_begin(); vit != part.vertices_end(); ++vit) {
        const double x = CGAL::to_double(vit->x());
        const double y = CGAL::to_double(vit->y());
        const auto it = idx_by_xy.find({x, y});
        if (it == idx_by_xy.end()) {
          throw std::runtime_error("CGAL partition produced a vertex not in input polygon");
        }
        vids.push_back(it->second);
      }

      for (int i = 1; i + 1 < m; ++i) {
        out_tris.push_back({vids[0], vids[i], vids[i + 1]});
      }
    }

    const double t1 = cpu_now_ms();
    const double elapsed_ms = (t1 - t0);

    if (static_cast<int>(out_tris.size()) != n - 2) {
      throw std::runtime_error("Unexpected triangle count (expected n-2)");
    }

    write_tri(args.output, pts, out_tris);
    std::cout << "hertel,vertices=" << pts.size()
              << ",triangles=" << out_tris.size()
              << ",time_ms=" << elapsed_ms << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }
}

