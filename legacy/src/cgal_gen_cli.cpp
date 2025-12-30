#include "polygon_io.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using Clock = std::chrono::high_resolution_clock;

namespace {

struct Args {
  std::string type;   // random_simple | convex
  int n = 0;
  std::uint64_t seed = 1;
  std::string output;
};

Args parse_args(int argc, char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    const std::string t = argv[i];
    if ((t == "--type") && i + 1 < argc) a.type = argv[++i];
    else if ((t == "--n") && i + 1 < argc) a.n = std::stoi(argv[++i]);
    else if ((t == "--seed") && i + 1 < argc) a.seed = static_cast<std::uint64_t>(std::stoull(argv[++i]));
    else if ((t == "-o" || t == "--output") && i + 1 < argc) a.output = argv[++i];
    else {
      throw std::runtime_error("Unknown or incomplete argument: " + t);
    }
  }
  if (a.type.empty()) throw std::runtime_error("Missing --type");
  if (a.n < 3) throw std::runtime_error("--n must be >= 3");
  if (a.output.empty()) throw std::runtime_error("Missing --output");
  return a;
}

void write_poly(const baselines::Polygon& poly, const std::string& path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Failed to open output file: " + path);
  out.setf(std::ios::fixed);
  out.precision(10);
  out << poly.size() << "\n";
  for (const auto& p : poly) out << p.x << " " << p.y << "\n";
}

}  // namespace

int main(int argc, char** argv) {
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Point = K::Point_2;

  try {
    const auto args = parse_args(argc, argv);

    std::mt19937_64 rng(args.seed);
    CGAL::Random cgal_rng(static_cast<unsigned int>(rng()));

    baselines::Polygon poly;
    poly.reserve(static_cast<std::size_t>(args.n));

    const auto t0 = Clock::now();

    if (args.type == "random_simple") {
      // Generate random points on a circle and polygonize with CGAL.
      // Using points on/near a circle avoids the pathological slow cases we observed
      // with points in a square for random_polygon_2 at larger n.
      std::vector<Point> pts;
      pts.reserve(static_cast<std::size_t>(args.n));
      // Points in a disc (not exactly on a circle) yield non-convex simple polygons in practice.
      CGAL::Random_points_in_disc_2<Point> gen(100.0, cgal_rng);
      for (int i = 0; i < args.n; ++i) pts.push_back(*gen++);

      std::vector<Point> out;
      out.reserve(static_cast<std::size_t>(args.n));
      CGAL::random_polygon_2(pts.size(), std::back_inserter(out), pts.begin());

      for (const auto& p : out) poly.push_back({CGAL::to_double(p.x()), CGAL::to_double(p.y())});
    } else if (args.type == "convex") {
      std::vector<Point> pts;
      pts.reserve(static_cast<std::size_t>(args.n));
      CGAL::Random_points_in_square_2<Point> gen(100.0, cgal_rng);
      for (int i = 0; i < args.n; ++i) pts.push_back(*gen++);

      std::vector<Point> hull;
      CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(hull));
      if (static_cast<int>(hull.size()) < 3) throw std::runtime_error("convex_hull_2 produced < 3 vertices");
      poly.reserve(hull.size());
      for (const auto& p : hull) poly.push_back({CGAL::to_double(p.x()), CGAL::to_double(p.y())});
    } else {
      throw std::runtime_error("Unknown --type: " + args.type);
    }

    const auto t1 = Clock::now();
    const double elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    write_poly(poly, args.output);
    std::cout << "cgal_gen,type=" << args.type << ",n=" << poly.size()
              << ",seed=" << args.seed << ",time_ms=" << elapsed_ms << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }
}


