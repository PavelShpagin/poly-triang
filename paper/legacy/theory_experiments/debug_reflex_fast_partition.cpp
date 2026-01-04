#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#define REFLEX_FAST_DEBUG
#include "../../methods/include/reflex_fast.hpp"

using reflex_fast::Point;

static std::vector<Point> make_radial_polygon(int n, uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> ang(0.0, 2.0 * M_PI);
  std::uniform_real_distribution<double> rad(50.0, 200.0);

  std::vector<double> a(n);
  for (int i = 0; i < n; ++i) a[i] = ang(rng);
  std::sort(a.begin(), a.end());

  std::vector<Point> pts(n);
  for (int i = 0; i < n; ++i) {
    const double r = rad(rng);
    pts[i].x = 500.0 + r * std::cos(a[i]);
    pts[i].y = 500.0 + r * std::sin(a[i]);
    pts[i].index = i;
  }
  return pts;
}

static double signed_area2(const std::vector<Point>& pts) {
  const int n = static_cast<int>(pts.size());
  double a = 0.0;
  for (int i = 0; i < n; ++i) {
    const int j = (i + 1) % n;
    a += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
  }
  return a;
}

struct LLV {
  int next = -1;
  int prev = -1;
  int orig = -1;
};

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

// PolyPartition-style split: duplicate both endpoints and rewire.
static void add_diagonal_split(std::vector<LLV>& ll, int& num, int i1, int i2) {
  const int new1 = num++;
  const int new2 = num++;

  ll[new1] = ll[i1];
  ll[new2] = ll[i2];

  // new2 follows i1, new1 follows i2
  ll[new2].next = ll[i2].next;
  ll[new1].next = ll[i1].next;

  ll[ll[i2].next].prev = new2;
  ll[ll[i1].next].prev = new1;

  ll[i1].next = new2;
  ll[new2].prev = i1;

  ll[i2].next = new1;
  ll[new1].prev = i2;
}

static int count_loops_after_splitting(int n, const std::vector<std::pair<int, int>>& diags) {
  const int maxv = n + 2 * static_cast<int>(diags.size());
  std::vector<LLV> ll(maxv);
  int num = n;

  for (int i = 0; i < n; ++i) {
    ll[i].next = (i + 1) % n;
    ll[i].prev = (i + n - 1) % n;
    ll[i].orig = i;
  }

  for (auto [a, b] : diags) {
    add_diagonal_split(ll, num, a, b);
  }

  std::vector<uint8_t> used(num, 0);
  int loops = 0;
  for (int s = 0; s < num; ++s) {
    if (used[s]) continue;
    ++loops;
    int cur = s;
    for (int safe = 0; safe < num + 10; ++safe) {
      used[cur] = 1;
      cur = ll[cur].next;
      if (cur == s) break;
    }
  }
  return loops;
}

static bool diagonals_look_planar(const std::vector<Point>& pts,
                                 const std::vector<std::pair<int, int>>& diags,
                                 int& bad_d0,
                                 int& bad_d1,
                                 int& bad_e0,
                                 int& bad_e1) {
  const int n = static_cast<int>(pts.size());

  // Diagonal-diagonal proper intersections (excluding shared endpoints).
  for (size_t i = 0; i < diags.size(); ++i) {
    const auto [a0, a1] = diags[i];
    for (size_t j = 0; j < i; ++j) {
      const auto [b0, b1] = diags[j];
      if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1) continue;
      if (proper_intersect(pts[a0], pts[a1], pts[b0], pts[b1])) {
        std::cout << "FOUND diagonal-diagonal intersection: "
                  << a0 << "-" << a1 << " with " << b0 << "-" << b1 << "\n";
        bad_d0 = a0;
        bad_d1 = a1;
        bad_e0 = b0;
        bad_e1 = b1;
        return false;
      }
    }
  }

  // Diagonal-boundary proper intersections (excluding incident edges).
  for (const auto [a0, a1] : diags) {
    for (int e = 0; e < n; ++e) {
      const int b0 = e;
      const int b1 = (e + 1) % n;
      if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1) continue;
      if (proper_intersect(pts[a0], pts[a1], pts[b0], pts[b1])) {
        std::cout << "FOUND diagonal-boundary intersection: diag "
                  << a0 << "-" << a1 << " with edge " << b0 << "-" << b1 << "\n";
        bad_d0 = a0;
        bad_d1 = a1;
        bad_e0 = b0;
        bad_e1 = b1;
        return false;
      }
    }
  }

  return true;
}

static double chain_x_at_y(const std::vector<Point>& pts,
                           const std::vector<int>& chain_verts,
                           double y) {
  const int len = static_cast<int>(chain_verts.size());
  int pos = 0;
  while (pos + 1 < len && pts[chain_verts[pos + 1]].y > y) ++pos;
  if (pos >= len - 1) pos = len - 2;
  if (pos < 0) pos = 0;
  const Point& p = pts[chain_verts[pos]];
  const Point& q = pts[chain_verts[pos + 1]];
  const double dy = p.y - q.y;
  if (std::abs(dy) < 1e-18) return std::min(p.x, q.x);
  return p.x + (p.y - y) / dy * (q.x - p.x);
}

static double seg_x_at_y(const Point& p, const Point& q, double y) {
  const double dy = p.y - q.y;
  if (std::abs(dy) < 1e-18) return std::min(p.x, q.x);
  const double t = (p.y - y) / dy;
  return p.x + t * (q.x - p.x);
}

int main(int argc, char** argv) {
  int n = 50;
  int seed = 0;
  if (argc >= 2) n = std::atoi(argv[1]);
  if (argc >= 3) seed = std::atoi(argv[2]);

  auto pts = make_radial_polygon(n, static_cast<uint32_t>(seed));
  std::cout << "input_area2=" << signed_area2(pts) << "\n";

  reflex_fast::FastTriangulator tri;
  tri.triangulate(pts);

  const auto& D = tri.diagonals();
  int bad_d0 = -1, bad_d1 = -1, bad_e0 = -1, bad_e1 = -1;
  const bool planar = diagonals_look_planar(pts, D, bad_d0, bad_d1, bad_e0, bad_e1);
  const int loops = count_loops_after_splitting(n, D);

  std::cout << "n=" << n << " seed=" << seed
            << " reflex=" << tri.reflex_count
            << " chains=" << tri.chain_count()
            << " diags=" << D.size()
            << " loops=" << loops
            << " expected_loops=" << (static_cast<int>(D.size()) + 1)
            << " tris=" << tri.triangles.size()
            << " expected_tris=" << (n - 2)
            << "\n";

  if (!planar) {
    std::cout << "BAD_DIAGONALS\n";
#ifdef REFLEX_FAST_DEBUG
    const int a = std::min(bad_d0, bad_d1);
    const int b = std::max(bad_d0, bad_d1);
    for (const auto& r : tri.debug_diag_records()) {
      if (r.a == a && r.b == b) {
        std::cout << "  record: reason=" << (r.reason ? r.reason : "<null>")
                  << " event_v=" << r.event_v
                  << " event_type=" << r.event_type
                  << " chosen_chain=" << r.chosen_chain
                  << " target_v=" << r.target_v
                  << "\n";
        // If this was a pending-based diagonal, print where pending was set.
        if (r.reason && std::string(r.reason).find("pending") != std::string::npos && r.chosen_chain >= 0) {
          for (const auto& pr : tri.debug_pending_records()) {
            if (pr.chain_id == r.chosen_chain && pr.pending_v == r.target_v) {
              std::cout << "  pending origin: chain_id=" << pr.chain_id
                        << " pending_v=" << pr.pending_v
                        << " set_at_event_v=" << pr.event_v
                        << " reason=" << (pr.reason ? pr.reason : "<null>")
                        << "\n";

              // Reconstruct the sweep state at the merge event that set this pending.
              const int mv = pr.event_v;
              double my = pts[mv].y + 1e-9;  // MERGE uses y+eps in our sweep.
              struct Item2 { double x; int cid; };
              std::vector<Item2> items2;
              for (int cid = 0; cid < tri.chain_count(); ++cid) {
                const auto ci = tri.debug_chain_info(cid);
                if (!(pts[ci.upper].y > my && my > pts[ci.lower].y)) continue;
                const auto cv = tri.debug_chain_vertices(cid);
                items2.push_back({chain_x_at_y(pts, cv, my), cid});
              }
              std::sort(items2.begin(), items2.end(), [](const Item2& A, const Item2& B) {
                if (A.x < B.x) return true;
                if (B.x < A.x) return false;
                return A.cid < B.cid;
              });
              std::cout << "  pending-set sweep_y=" << my << " at merge v=" << mv
                        << " x=" << pts[mv].x << "\n";
              std::cout << "  active chains at pending-set sweep_y (x,cid):\n";
              for (const auto& it2 : items2) {
                std::cout << "    x=" << it2.x << " cid=" << it2.cid << "\n";
              }

              break;
            }
          }
        }
        if (r.chosen_chain >= 0) {
          const auto info = tri.debug_chain_info(r.chosen_chain);
          const auto verts = tri.debug_chain_vertices(r.chosen_chain);
          std::cout << "  chain[" << info.id << "]: upper=" << info.upper
                    << " lower=" << info.lower
                    << " pending=" << info.pending
                    << " cached_pos=" << info.cached_pos
                    << " len=" << verts.size() << "\n";
          if (info.cached_pos >= 0 && info.cached_pos + 1 < static_cast<int>(verts.size())) {
            const int u = verts[info.cached_pos];
            const int w = verts[info.cached_pos + 1];
            std::cout << "    curr edge: (" << u << ")->(" << w << ")\n";
            std::cout << "      u=(" << pts[u].x << "," << pts[u].y << ")\n";
            std::cout << "      w=(" << pts[w].x << "," << pts[w].y << ")\n";
          }
          std::cout << "  event v=(" << pts[r.event_v].x << "," << pts[r.event_v].y << ")\n";
          std::cout << "  target t=(" << pts[r.target_v].x << "," << pts[r.target_v].y << ")\n";
          const int lo = std::max(0, std::min({r.event_v, r.target_v, bad_e0, bad_e1}) - 3);
          const int hi = std::min(n - 1, std::max({r.event_v, r.target_v, bad_e0, bad_e1}) + 3);
          std::cout << "  local vertices [" << lo << ".." << hi << "] (cyclic order indices):\n";
          for (int i = lo; i <= hi; ++i) {
            std::cout << "    " << i << ": (" << pts[i].x << "," << pts[i].y << ")\n";
          }

          // Dump active-ish chain ordering at the event y (best-effort reconstruction).
          double y = pts[r.event_v].y;
          if (r.event_type == 0 || r.event_type == 2) y -= 1e-9;  // START/SPLIT
          else y += 1e-9;                                         // END/MERGE
          const double qx = pts[r.event_v].x;
          struct Item { double x; int cid; };
          std::vector<Item> items;
          for (int cid = 0; cid < tri.chain_count(); ++cid) {
            const auto ci = tri.debug_chain_info(cid);
            // Active if y lies between upper and lower.
            if (!(pts[ci.upper].y > y && y > pts[ci.lower].y)) continue;
            const auto cv = tri.debug_chain_vertices(cid);
            const double cx = chain_x_at_y(pts, cv, y);
            items.push_back({cx, cid});
          }
          std::sort(items.begin(), items.end(), [](const Item& A, const Item& B) {
            if (A.x < B.x) return true;
            if (B.x < A.x) return false;
            return A.cid < B.cid;
          });
          std::cout << "  sweep_y=" << y << " query_x=" << qx << "\n";
          std::cout << "  active chains at sweep_y (x,cid):\n";
          for (const auto& it : items) {
            std::cout << "    x=" << it.x << " cid=" << it.cid << "\n";
          }

          // Compare with the true left boundary edge among ALL polygon edges at this y.
          struct EdgeX { double x; int i; int j; };
          std::vector<EdgeX> ex;
          ex.reserve(n);
          for (int i = 0; i < n; ++i) {
            const int j = (i + 1) % n;
            const double yi = pts[i].y;
            const double yj = pts[j].y;
            // Strict crossing of horizontal line (ignore if endpoint lies on the line).
            if ((yi > y && yj < y) || (yj > y && yi < y)) {
              ex.push_back({seg_x_at_y(pts[i], pts[j], y), i, j});
            }
          }
          std::sort(ex.begin(), ex.end(), [](const EdgeX& A, const EdgeX& B) { return A.x < B.x; });
          int best = -1;
          for (int k = 0; k < static_cast<int>(ex.size()); ++k) {
            if (ex[k].x < qx) best = k;
            else break;
          }
          std::cout << "  boundary crossings at sweep_y: " << ex.size() << "\n";
          if (best >= 0) {
            const auto& e = ex[best];
            std::cout << "  true left edge: (" << e.i << ")->(" << e.j << ") x=" << e.x << "\n";
            std::cout << "    pi=(" << pts[e.i].x << "," << pts[e.i].y << ")\n";
            std::cout << "    pj=(" << pts[e.j].x << "," << pts[e.j].y << ")\n";
          } else {
            std::cout << "  true left edge: <none>\n";
          }
        }
        break;
      }
    }
#endif
    return 1;
  }
  if (loops != static_cast<int>(D.size()) + 1) {
    std::cout << "Diagonals:\n";
    for (auto [a, b] : D) {
      std::cout << "  " << a << " " << b << "\n";
    }
    std::cout << "BAD_PARTITION\n";
    return 1;
  }
  if (static_cast<int>(tri.triangles.size()) != n - 2) {
    std::cout << "BAD_TRI_COUNT\n";
    return 2;
  }
  std::cout << "OK\n";
  return 0;
}


