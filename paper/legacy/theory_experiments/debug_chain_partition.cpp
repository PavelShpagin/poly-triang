#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../../methods/include/reflex_chain_triangulate.hpp"

using reflex_tri::Point;

static int orient(const Point& a, const Point& b, const Point& c) {
  const long double bax = static_cast<long double>(b.x) - static_cast<long double>(a.x);
  const long double bay = static_cast<long double>(b.y) - static_cast<long double>(a.y);
  const long double cax = static_cast<long double>(c.x) - static_cast<long double>(a.x);
  const long double cay = static_cast<long double>(c.y) - static_cast<long double>(a.y);
  const long double v = bax * cay - bay * cax;
  if (v > 0) return 1;
  if (v < 0) return -1;
  return 0;
}

static bool proper_intersect(const Point& a, const Point& b, const Point& c, const Point& d) {
  const int o1 = orient(a, b, c);
  const int o2 = orient(a, b, d);
  const int o3 = orient(c, d, a);
  const int o4 = orient(c, d, b);
  return (o1 * o2 < 0) && (o3 * o4 < 0);
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <polygon.poly>\n";
    return 1;
  }

  const std::string path = argv[1];
  std::ifstream fin(path);
  if (!fin) {
    std::cerr << "Error: cannot open " << path << "\n";
    return 1;
  }

  int n = 0;
  fin >> n;
  if (n < 3) {
    std::cerr << "Error: n<3\n";
    return 1;
  }

  std::vector<Point> pts(n);
  for (int i = 0; i < n; ++i) {
    fin >> pts[i].x >> pts[i].y;
    pts[i].index = i;
  }

  reflex_tri::Triangulator tri;
  std::string err;
  try {
    (void)tri.triangulate(pts);
  } catch (const std::exception& e) {
    err = e.what();
  }

  const auto& diags = tri.debug_diagonals();
  std::cerr << "n=" << n << " diagonals=" << diags.size() << " reflex=" << tri.reflex_count() << "\n";
  if (!err.empty()) {
    std::cerr << "triangulate() threw: " << err << "\n";
  }

  // Check diagonal-boundary intersections.
  for (size_t i = 0; i < diags.size(); ++i) {
    const int a = diags[i].first;
    const int b = diags[i].second;
    const Point& pa = pts[a];
    const Point& pb = pts[b];
    for (int e = 0; e < n; ++e) {
      const int c = e;
      const int d = (e + 1) % n;
      if (c == a || c == b || d == a || d == b) continue;
      if (proper_intersect(pa, pb, pts[c], pts[d])) {
        std::cerr << "DIAG-BOUNDARY INTERSECTION: diag(" << a << "," << b << ") x edge("
                  << c << "," << d << ")\n";
        for (const auto& r : tri.debug_diag_records()) {
          if (r.a == a && r.b == b) {
            std::cerr << "  record: event_v=" << r.event_v
                      << " event_type=" << r.event_type
                      << " chosen_chain=" << r.chosen_chain
                      << " target_v=" << r.target_v
                      << " reason=" << (r.reason ? r.reason : "") << "\n";
            // Also show where the pending target was set.
            for (const auto& pr : tri.debug_pending_records()) {
              if (pr.pending_v == r.target_v) {
                std::cerr << "  pending-set: chain_id=" << pr.chain_id
                          << " pending_v=" << pr.pending_v
                          << " event_v=" << pr.event_v
                          << " reason=" << (pr.reason ? pr.reason : "") << "\n";
              }
            }
            break;
          }
        }
        return 2;
      }
    }
  }

  // Check diagonal-diagonal intersections.
  for (size_t i = 0; i < diags.size(); ++i) {
    const int a = diags[i].first;
    const int b = diags[i].second;
    for (size_t j = i + 1; j < diags.size(); ++j) {
      const int c = diags[j].first;
      const int d = diags[j].second;
      // skip shared endpoints
      if (a == c || a == d || b == c || b == d) continue;
      if (proper_intersect(pts[a], pts[b], pts[c], pts[d])) {
        std::cerr << "DIAG-DIAG INTERSECTION: diag(" << a << "," << b << ") x diag(" << c << "," << d << ")\n";
        return 3;
      }
    }
  }

  std::cerr << "No proper intersections detected.\n";
  return 0;
}


