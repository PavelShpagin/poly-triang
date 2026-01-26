#pragma once
/**
 * Chain-based polygon triangulation (paper implementation).
 *
 * Target complexity: O(n + k log k), where k is the number of local maxima
 * (equivalently, local minima) with respect to the sweep direction.
 *
 * Notes:
 * - We assume a simple polygon without holes.
 * - We support equal y-coordinates via a deterministic tie-break (y, then x),
 *   but the theoretical paper assumes general position.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <stdexcept>
#include <ostream>
#include <sstream>
#include <vector>

// Optional boundary-crossing validation for proposed diagonals.
//
// When enabled, we reject diagonals that intersect the polygon boundary (excluding
// incident edges). This is a defensive guard for near-degenerate inputs.
//
// NOTE: This is off by default for performance; it can dominate runtime on high-k
// polygons.
//
// Build-time toggle:
// - 0 (default): disabled (fast)
// - 1: enabled (defensive)
#ifndef REFLEX_ENABLE_DIAGONAL_CROSS_CHECK
#define REFLEX_ENABLE_DIAGONAL_CROSS_CHECK 0
#endif

namespace reflex_tri {

struct Point {
  double x = 0.0;
  double y = 0.0;
  int index = -1;
};

struct Triangle {
  int v0, v1, v2;
};

namespace detail {

static constexpr double kEpsGeom = 1e-10;
static constexpr double kEpsMath = 1e-16;

inline double nextafter_up(double x) {
  // Next representable double toward +infinity (exact, no libm call).
  if (std::isnan(x) || x == std::numeric_limits<double>::infinity()) return x;
  if (x == 0.0) return std::numeric_limits<double>::denorm_min();
  std::uint64_t u = 0;
  std::memcpy(&u, &x, sizeof(u));
  if (x > 0.0) {
    ++u;
  } else {
    --u;
  }
  double out = 0.0;
  std::memcpy(&out, &u, sizeof(out));
  return out;
}

inline double nextafter_down(double x) {
  // Next representable double toward -infinity (exact, no libm call).
  if (std::isnan(x) || x == -std::numeric_limits<double>::infinity()) return x;
  if (x == 0.0) return -std::numeric_limits<double>::denorm_min();
  std::uint64_t u = 0;
  std::memcpy(&u, &x, sizeof(u));
  if (x > 0.0) {
    --u;
  } else {
    ++u;
  }
  double out = 0.0;
  std::memcpy(&out, &u, sizeof(out));
  return out;
}

inline double cross(const Point& a, const Point& b, const Point& c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

inline bool is_convex_polygon(const std::vector<Point>& pts) {
  const int n = static_cast<int>(pts.size());
  if (n < 4) return true;
  int sign = 0;  // 0 = unknown, +1 = left turns, -1 = right turns
  for (int i = 0; i < n; ++i) {
    const int p = (i == 0) ? (n - 1) : (i - 1);
    const int nx = (i + 1 == n) ? 0 : (i + 1);
    const double c = cross(pts[p], pts[i], pts[nx]);
    if (c > kEpsGeom) {
      if (sign < 0) return false;
      sign = 1;
    } else if (c < -kEpsGeom) {
      if (sign > 0) return false;
      sign = -1;
    }
  }
  return true;
}

inline double signed_area(const std::vector<Point>& pts) {
  double a = 0.0;
  const int n = static_cast<int>(pts.size());
  for (int i = 0; i < n; ++i) {
    const int j = (i + 1) % n;
    a += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
  }
  return 0.5 * a;
}

inline bool below(const std::vector<Point>& pts, int a, int b) {
  // Strict, deterministic order:
  // - primary: y (smaller y is below)
  // - tie: x (larger x is below)  [matches PolyPartition convention]
  // - final tie: index (larger index is below) to guarantee strict ordering
  if (std::abs(pts[a].y - pts[b].y) > kEpsGeom) {
    return pts[a].y < pts[b].y;
  }
  if (pts[a].x > pts[b].x) return true;
  if (pts[a].x < pts[b].x) return false;
  return a > b;
}

inline bool above(const std::vector<Point>& pts, int a, int b) {
  return below(pts, b, a);
}

inline bool is_reflex_vertex_ccw(const std::vector<Point>& pts, int i) {
  const int n = static_cast<int>(pts.size());
  const int p = (i - 1 + n) % n;
  const int nx = (i + 1) % n;
  return cross(pts[p], pts[i], pts[nx]) < -kEpsGeom;
}

}  // namespace detail

class Triangulator {
public:
  const std::vector<Triangle>& triangulate(std::vector<Point>& pts) {
    triangles_.clear();
    diagonals_.clear();
    reflex_count_ = 0;
    dbg_faces_ = 0;
    dbg_failed_walks_ = 0;
    dbg_degenerate_faces_ = 0;
    dbg_iter_limit_faces_ = 0;
    dbg_cw_faces_ = 0;
    dbg_unused_halfedges_ = 0;
    dbg_euler_mismatch_ = 0;
    dbg_total_walks_ = 0;
    dbg_sum_face_verts_ = 0;
    dbg_half_edges_ = 0;
    dbg_used_half_edges_ = 0;

    const int n = static_cast<int>(pts.size());
    if (n < 3) return triangles_;
    // We always output exactly (n-2) triangles on success.
    triangles_.reserve(static_cast<std::size_t>(n - 2));

    // Note: we do not do a separate convexity pre-check here. Convex polygons are
    // detected via reflex_count_ == 0 after orientation normalization below.

    // Ensure CCW orientation (algorithm assumes CCW).
    if (detail::signed_area(pts) < 0.0) std::reverse(pts.begin(), pts.end());

    // NOTE: The paper assumes general position (distinct y-coordinates).
    // For performance (and to preserve the output-sensitive bound), we do not
    // enforce this via an O(n log n) check here.

    // Count all reflex vertices (paper's r) and cache convexity flags for reuse in
    // vertex classification (avoids recomputing cross products in Phase 1).
    if (static_cast<int>(is_convex_buf_.size()) != n) {
      is_convex_buf_.assign(n, 0);
    } else {
      std::fill(is_convex_buf_.begin(), is_convex_buf_.end(), 0);
    }
    for (int i = 0; i < n; ++i) {
      const int p = (i - 1 + n) % n;
      const int nx = (i + 1) % n;
      const double c = detail::cross(pts[p], pts[i], pts[nx]);
      if (c > detail::kEpsGeom) is_convex_buf_[i] = 1;
      if (c < -detail::kEpsGeom) ++reflex_count_;
    }
    // Roughly, the number of decomposition diagonals is O(reflex_count_).
    diagonals_.reserve(static_cast<std::size_t>(reflex_count_) + 8);

    if (n == 3) {
      triangles_.push_back({0, 1, 2});
      return triangles_;
    }

    if (reflex_count_ == 0) {
      // Convex -> fan triangulation.
      triangles_.resize(n - 2);
      for (int i = 0; i < n - 2; ++i) triangles_[i] = {0, i + 1, i + 2};
      return triangles_;
    }

#ifdef REFLEX_TRI_TIMING
    using clock = std::chrono::high_resolution_clock;
    const auto t0 = clock::now();
#endif

    // Build a spatial index over boundary edges for diagonal validation.
#if REFLEX_ENABLE_DIAGONAL_CROSS_CHECK
    build_boundary_bvh(pts);
#endif

    build_vertex_types_and_chains(pts);

#ifdef REFLEX_TRI_TIMING
    const auto t1 = clock::now();
#endif

    decompose_into_monotone(pts);

#ifdef REFLEX_TRI_TIMING
    const auto t2 = clock::now();
#endif

    triangulate_monotone_faces(pts);

#ifdef REFLEX_TRI_TIMING
    const auto t3 = clock::now();
    phase1_ms_ = std::chrono::duration<double, std::milli>(t1 - t0).count();
    phase2_ms_ = std::chrono::duration<double, std::milli>(t2 - t1).count();
    phase3_ms_ = std::chrono::duration<double, std::milli>(t3 - t2).count();
#endif

    // Sanity: simple polygon triangulation => n-2 triangles.
    // (If a baseline rejects a degenerate polygon, we fail rather than "fallback".)
    if (static_cast<int>(triangles_.size()) != n - 2) {
      // Do not attempt any fallback; surface the issue.
      throw std::runtime_error(
          "triangulation failed: triangle count != n-2 "
          "(n=" + std::to_string(n) +
          ", triangles=" + std::to_string(triangles_.size()) +
          ", diagonals=" + std::to_string(diagonals_.size()) +
          ", faces=" + std::to_string(dbg_faces_) +
          ", failed_walks=" + std::to_string(dbg_failed_walks_) +
          ", degenerate=" + std::to_string(dbg_degenerate_faces_) +
          ", iter_limit=" + std::to_string(dbg_iter_limit_faces_) +
          ", cw_faces=" + std::to_string(dbg_cw_faces_) +
          ", unused_he=" + std::to_string(dbg_unused_halfedges_) +
          ", euler_miss=" + std::to_string(dbg_euler_mismatch_) +
          ", total_walks=" + std::to_string(dbg_total_walks_) +
          ", sum_face_verts=" + std::to_string(dbg_sum_face_verts_) +
          ", half_edges=" + std::to_string(dbg_half_edges_) +
          ", used_half_edges=" + std::to_string(dbg_used_half_edges_) +
          ", implied_tri=" + std::to_string(
              (dbg_faces_ > 0) ? (dbg_sum_face_verts_ - 2LL * dbg_faces_) : 0LL) +
          ")");
    }

    return triangles_;
  }

  int reflex_count() const { return reflex_count_; }
  const std::vector<std::pair<int, int>>& debug_diagonals() const { return diagonals_; }
  const std::vector<Triangle>& debug_triangles() const { return triangles_; }
#ifdef REFLEX_TRI_TIMING
  double phase1_ms() const { return phase1_ms_; }
  double phase2_ms() const { return phase2_ms_; }
  double phase3_ms() const { return phase3_ms_; }
#endif
  struct DiagRecord {
    int a = -1;
    int b = -1;
    int event_v = -1;
    int event_type = -1;   // casted from VType
    int chosen_chain = -1;
    int target_v = -1;
    const char* reason = "";
  };
  const std::vector<DiagRecord>& debug_diag_records() const { return diag_records_; }
  struct PendingRecord {
    int chain_id = -1;
    int pending_v = -1;
    int event_v = -1;
    const char* reason = "";
  };
  const std::vector<PendingRecord>& debug_pending_records() const { return pending_records_; }

private:
  enum class VType { Start, End, Split, Merge, Regular };
  struct Event;

  // Fast open-addressing set for diagonal keys (uint64_t), to avoid per-insert heap
  // allocations of std::unordered_set in the critical path.
  struct DiagSet {
    std::vector<std::uint64_t> tab;
    std::size_t mask = 0;

    static std::uint64_t mix(std::uint64_t x) {
      // splitmix64 mix
      x ^= x >> 33;
      x *= 0xff51afd7ed558ccdULL;
      x ^= x >> 33;
      x *= 0xc4ceb9fe1a85ec53ULL;
      x ^= x >> 33;
      return x;
    }

    void reset(std::size_t expected) {
      std::size_t cap = 1;
      const std::size_t need = std::max<std::size_t>(16, expected * 2);
      while (cap < need) cap <<= 1;
      tab.assign(cap, 0);
      mask = cap - 1;
    }

    bool insert(std::uint64_t key) {
      // 0 is the empty sentinel; diag_key never returns 0 for u!=v.
      std::size_t i = static_cast<std::size_t>(mix(key)) & mask;
      while (true) {
        std::uint64_t& slot = tab[i];
        if (slot == 0) {
          slot = key;
          return true;
        }
        if (slot == key) return false;
        i = (i + 1) & mask;
      }
    }
  };

  struct Chain {
    std::vector<int> verts;  // from upper (local max) to lower (local min)
    int curr = 0;            // current edge is verts[curr] -> verts[curr+1]
    // Pending merge vertex (paper algorithm). If non-null, this merge vertex
    // should be connected downward at the next suitable extremum.
    int pending = -1;
    int upper = -1;
    int lower = -1;
    double cache_x = 0.0;
    int cache_stamp = -1;
    // Current edge parameters for fast x_at(y):
    // If edge is non-horizontal, x(y) = base_x - slope * y (derived from the edge endpoints).
    int edge_a = -1;
    int edge_b = -1;
    double edge_base_x = 0.0;
    double edge_slope = 0.0;

    void reset_runtime() {
      curr = 0;
      pending = -1;
      cache_x = 0.0;
      cache_stamp = -1;
      edge_a = -1;
      edge_b = -1;
      edge_base_x = 0.0;
      edge_slope = 0.0;
    }

    void advance(double y, const std::vector<Point>& pts) {
      // Move down the chain until current edge spans y.
      const int old_curr = curr;
      while (curr + 1 < static_cast<int>(verts.size()) &&
             pts[verts[curr + 1]].y > y + detail::kEpsGeom) {
        ++curr;
      }
      if (curr + 1 >= static_cast<int>(verts.size())) {
        curr = static_cast<int>(verts.size()) - 2;
      }
      if (curr < 0) curr = 0;
      if (curr != old_curr) {
        cache_stamp = -1;
        edge_a = -1;
      }
    }

    int slab_entry_vertex() const {
      // Upper vertex of current edge.
      return verts[curr];
    }

    // NOTE: x-at-sweep is computed by the owning Triangulator so it can use
    // the shared sweep_y_ cache and precomputed edge parameters.
  };

  // Treap keyed by chain order at current sweep_y (see paper).
  struct TreapNode {
    int chain_id = -1;
    std::uint32_t prio = 0;
    TreapNode* left = nullptr;
    TreapNode* right = nullptr;
  };

  std::vector<VType> types_;
  std::vector<Chain> chains_;
  std::vector<int> max_to_chain_;  // [vertex] -> left-boundary chain_id (for local maxima)
  std::vector<int> min_to_chain_;  // [vertex] -> left-boundary chain_id (for local minima)

  double sweep_y_ = 0.0;
  double sweep_y_eps_ = 0.0;
  int sweep_stamp_ = 0;
  std::vector<std::pair<int, int>> diagonals_;
  // For uniqueness without an O(|D| log |D|) sort pass.
  DiagSet diag_set_;
  std::vector<DiagRecord> diag_records_;
  std::vector<PendingRecord> pending_records_;
  std::vector<Triangle> triangles_;
  int reflex_count_ = 0;
  std::vector<unsigned char> is_convex_buf_;
  std::vector<Event> events_;
  // Phase 3 scratch (kept as members to avoid per-call heap deallocations).
  std::vector<int> adj_deg_;
  std::vector<int> adj_off_;
  std::vector<int> adj_flat_;
  std::vector<int> adj_flat_tmp_;
  std::vector<int> adj_cur_;
  std::vector<int> halfedge_src_;
  std::vector<int> halfedge_key_;
  std::vector<int> halfedge_idx1_;
  std::vector<int> halfedge_idx2_;
  std::vector<unsigned char> halfedge_used_;
  std::vector<int> halfedge_rev_e_;
  std::vector<int> count_buf_;
  std::vector<int> face_buf_;
#ifdef REFLEX_TRI_TIMING
  double phase1_ms_ = 0.0;
  double phase2_ms_ = 0.0;
  double phase3_ms_ = 0.0;
#endif
  // Debug stats from face extraction (helpful to diagnose degeneracy / bugs).
  int dbg_faces_ = 0;
  int dbg_failed_walks_ = 0;
  int dbg_degenerate_faces_ = 0;
  int dbg_iter_limit_faces_ = 0;
  int dbg_cw_faces_ = 0;
  int dbg_unused_halfedges_ = 0;
  int dbg_euler_mismatch_ = 0;
  int dbg_total_walks_ = 0;
  
  // Temporary reference to pts for diagonal validation in add_diagonal_unique.
  // Kept as a pointer to avoid O(n) copies.
  const std::vector<Point>* diag_check_pts_ = nullptr;

  // Boundary-edge BVH for fast "diagonal crosses boundary?" queries.
  struct Aabb {
    double minx = 0.0;
    double miny = 0.0;
    double maxx = 0.0;
    double maxy = 0.0;
  };
  struct BvhNode {
    Aabb box;
    int left = -1;
    int right = -1;
    int start = 0;  // index into boundary_edge_ids_
    int count = 0;  // leaf size; 0 for internal nodes
  };
  std::vector<int> boundary_edge_ids_;   // edge id e corresponds to boundary edge (e, (e+1)%n)
  std::vector<BvhNode> boundary_bvh_;
  int boundary_bvh_root_ = -1;
  long long dbg_sum_face_verts_ = 0;
  std::size_t dbg_half_edges_ = 0;
  std::size_t dbg_used_half_edges_ = 0;

  // Scratch buffers (reused across faces) to avoid per-face allocations in
  // monotone triangulation.
  std::vector<char> is_left_buf_;
  std::vector<int> is_left_touched_;
  std::vector<int> tmp_chain1_;
  std::vector<int> tmp_chain2_;
  std::vector<int> tmp_sorted_;
  std::vector<int> tmp_stack_;
  std::vector<int> ear_clip_list_;

  // --- Treap utilities ---
  std::uint32_t rng_ = 0xA341316Cu;
  std::uint32_t next_prio() {
    // xorshift32
    std::uint32_t x = rng_;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng_ = x;
    return x;
  }

  // =========================
  // Status structure:
  // randomized treap as a TRUE BST in the left-to-right order of active chains
  // (ordered by x_at(sweep_y_) with tie-breaker by chain_id).
  // =========================
  struct StatusNode {
    int chain_id = -1;
    std::uint32_t prio = 0;
    StatusNode* left = nullptr;
    StatusNode* right = nullptr;
  };

  std::vector<StatusNode*> status_node_of_chain_;
  std::vector<StatusNode> status_pool_;

  static std::uint64_t diag_key(int a, int b) {
    // Undirected key with a<b.
    if (a > b) std::swap(a, b);
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) |
           static_cast<std::uint32_t>(b);
  }

  static Aabb aabb_from_points(const Point& a, const Point& b) {
    Aabb bb;
    bb.minx = std::min(a.x, b.x);
    bb.maxx = std::max(a.x, b.x);
    bb.miny = std::min(a.y, b.y);
    bb.maxy = std::max(a.y, b.y);
    return bb;
  }

  static bool aabb_overlaps(const Aabb& a, const Aabb& b) {
    if (a.maxx < b.minx || b.maxx < a.minx) return false;
    if (a.maxy < b.miny || b.maxy < a.miny) return false;
    return true;
  }

  static bool segment_intersects_aabb(const Point& p0, const Point& p1, const Aabb& box) {
    // Slab test on [0,1] segment parameter.
    double t0 = 0.0;
    double t1 = 1.0;

    const double dx = p1.x - p0.x;
    if (std::abs(dx) < detail::kEpsMath) {
      if (p0.x < box.minx || p0.x > box.maxx) return false;
    } else {
      const double inv = 1.0 / dx;
      double ta = (box.minx - p0.x) * inv;
      double tb = (box.maxx - p0.x) * inv;
      if (ta > tb) std::swap(ta, tb);
      t0 = std::max(t0, ta);
      t1 = std::min(t1, tb);
      if (t0 > t1) return false;
    }

    const double dy = p1.y - p0.y;
    if (std::abs(dy) < detail::kEpsMath) {
      if (p0.y < box.miny || p0.y > box.maxy) return false;
    } else {
      const double inv = 1.0 / dy;
      double ta = (box.miny - p0.y) * inv;
      double tb = (box.maxy - p0.y) * inv;
      if (ta > tb) std::swap(ta, tb);
      t0 = std::max(t0, ta);
      t1 = std::min(t1, tb);
      if (t0 > t1) return false;
    }

    return true;
  }

  Aabb edge_aabb(const std::vector<Point>& pts, int e) const {
    const int n = static_cast<int>(pts.size());
    const int a = e;
    const int b = (e + 1 == n) ? 0 : (e + 1);
    return aabb_from_points(pts[a], pts[b]);
  }

  double edge_mid_coord(const std::vector<Point>& pts, int e, int axis) const {
    const int n = static_cast<int>(pts.size());
    const int a = e;
    const int b = (e + 1 == n) ? 0 : (e + 1);
    if (axis == 0) return 0.5 * (pts[a].x + pts[b].x);
    return 0.5 * (pts[a].y + pts[b].y);
  }

  static Aabb aabb_union(const Aabb& a, const Aabb& b) {
    Aabb out;
    out.minx = std::min(a.minx, b.minx);
    out.miny = std::min(a.miny, b.miny);
    out.maxx = std::max(a.maxx, b.maxx);
    out.maxy = std::max(a.maxy, b.maxy);
    return out;
  }

  int build_boundary_bvh_rec(const std::vector<Point>& pts, int start, int end) {
    const int node_idx = static_cast<int>(boundary_bvh_.size());
    boundary_bvh_.push_back(BvhNode{});

    // Compute node bbox.
    Aabb bb;
    bool bb_init = false;
    for (int i = start; i < end; ++i) {
      const Aabb ebox = edge_aabb(pts, boundary_edge_ids_[i]);
      if (!bb_init) {
        bb = ebox;
        bb_init = true;
      } else {
        bb = aabb_union(bb, ebox);
      }
    }
    boundary_bvh_[node_idx].box = bb;

    const int count = end - start;
    static constexpr int kLeafSize = 16;
    if (count <= kLeafSize) {
      boundary_bvh_[node_idx].start = start;
      boundary_bvh_[node_idx].count = count;
      return node_idx;
    }

    // Split by the widest axis.
    const double dx = bb.maxx - bb.minx;
    const double dy = bb.maxy - bb.miny;
    const int axis = (dx >= dy) ? 0 : 1;
    const int mid = start + count / 2;

    auto it0 = boundary_edge_ids_.begin() + start;
    auto itM = boundary_edge_ids_.begin() + mid;
    auto it1 = boundary_edge_ids_.begin() + end;
    std::nth_element(it0, itM, it1, [&](int ea, int eb) {
      const double ca = edge_mid_coord(pts, ea, axis);
      const double cb = edge_mid_coord(pts, eb, axis);
      if (ca < cb) return true;
      if (cb < ca) return false;
      return ea < eb;
    });

    const int left = build_boundary_bvh_rec(pts, start, mid);
    const int right = build_boundary_bvh_rec(pts, mid, end);
    boundary_bvh_[node_idx].left = left;
    boundary_bvh_[node_idx].right = right;
    boundary_bvh_[node_idx].count = 0;
    return node_idx;
  }

  void build_boundary_bvh(const std::vector<Point>& pts) {
#if REFLEX_ENABLE_DIAGONAL_CROSS_CHECK
    const int n = static_cast<int>(pts.size());
    boundary_edge_ids_.resize(n);
    for (int i = 0; i < n; ++i) boundary_edge_ids_[i] = i;
    boundary_bvh_.clear();
    boundary_bvh_.reserve(static_cast<std::size_t>(2 * n + 1));
    boundary_bvh_root_ = build_boundary_bvh_rec(pts, 0, n);
#else
    (void)pts;
    boundary_edge_ids_.clear();
    boundary_bvh_.clear();
    boundary_bvh_root_ = -1;
#endif
  }

  static bool proper_segment_cross(const Point& a, const Point& b, const Point& c, const Point& d) {
    // Proper intersection test (excludes collinear/touching cases).
    const double o1 = detail::cross(a, b, c);
    const double o2 = detail::cross(a, b, d);
    const double o3 = detail::cross(c, d, a);
    const double o4 = detail::cross(c, d, b);
    const double eps = detail::kEpsGeom;
    if ((o1 > eps && o2 < -eps) || (o1 < -eps && o2 > eps)) {
      if ((o3 > eps && o4 < -eps) || (o3 < -eps && o4 > eps)) return true;
    }
    return false;
  }

  bool diagonal_crosses_boundary_bvh(const std::vector<Point>& pts, int node_idx, const Point& a, const Point& b, int u, int v) const {
    const BvhNode& node = boundary_bvh_[node_idx];
    if (!segment_intersects_aabb(a, b, node.box)) return false;
    if (node.left < 0 && node.right < 0) {
      const int n = static_cast<int>(pts.size());
      for (int i = 0; i < node.count; ++i) {
        const int e = boundary_edge_ids_[node.start + i];
        const int p = e;
        const int q = (e + 1 == n) ? 0 : (e + 1);
        // Skip edges incident to u or v.
        if (p == u || p == v || q == u || q == v) continue;
        if (proper_segment_cross(a, b, pts[p], pts[q])) return true;
      }
      return false;
    }
    if (node.left >= 0 && diagonal_crosses_boundary_bvh(pts, node.left, a, b, u, v)) return true;
    if (node.right >= 0 && diagonal_crosses_boundary_bvh(pts, node.right, a, b, u, v)) return true;
    return false;
  }

  // Check if segment (u,v) crosses any polygon boundary edge.
  // This ensures the diagonal doesn't exit the polygon.
  bool diagonal_crosses_boundary(const std::vector<Point>& pts, int u, int v) const {
    if (boundary_bvh_root_ < 0) return false;
    const Point& a = pts[u];
    const Point& b = pts[v];
    return diagonal_crosses_boundary_bvh(pts, boundary_bvh_root_, a, b, u, v);
  }

  bool add_diagonal_unique(int u, int v, int event_v, VType event_type, int chosen_chain,
                           int target_v, const char* reason) {
    const int n = static_cast<int>(types_.size());
    // Never insert a boundary edge as a "diagonal" (can create duplicate edges in Phase 3).
    if (n > 0) {
      if ((u + 1) % n == v || (v + 1) % n == u) return false;
    }
    if (u == v) return false;
    if (u > v) std::swap(u, v);
    
    // Optional defensive check (O(n) scan). Keep off for benchmarks.
#if REFLEX_ENABLE_DIAGONAL_CROSS_CHECK
    if (diag_check_pts_ && diagonal_crosses_boundary(*diag_check_pts_, u, v)) return false;
#endif
    const std::uint64_t key = diag_key(u, v);
    if (!diag_set_.insert(key)) return false;
    diagonals_.push_back({u, v});
#ifndef NDEBUG
    DiagRecord r;
    r.a = u;
    r.b = v;
    r.event_v = event_v;
    r.event_type = static_cast<int>(event_type);
    r.chosen_chain = chosen_chain;
    r.target_v = target_v;
    r.reason = reason;
    diag_records_.push_back(r);
#else
    (void)event_v;
    (void)event_type;
    (void)chosen_chain;
    (void)target_v;
    (void)reason;
#endif
    return true;
  }

  void advance_chain(int cid, const std::vector<Point>& pts) {
    Chain& ch = chains_[cid];
    const int n = static_cast<int>(pts.size());
    const int old_curr = ch.curr;
    // Step the chain pointer down while the next vertex is still above sweep_y_.
    while (ch.curr + 1 < static_cast<int>(ch.verts.size()) &&
           pts[ch.verts[ch.curr + 1]].y > sweep_y_eps_) {
      ++ch.curr;
      const int vreg = ch.verts[ch.curr];
      if (types_[vreg] != VType::Regular) continue;
      // REGULAR/left-chain case (PolyPartition): Below(v_i, v_{i-1}).
      // These vertices can resolve a pending MERGE helper in O(1) amortized time
      // without being explicit sweep events.
      const int prev = (vreg - 1 + n) % n;
      if (!detail::below(pts, vreg, prev)) continue;

      if (ch.pending != -1) {
        add_diagonal_unique(vreg, ch.pending, /*event_v=*/vreg, VType::Regular,
                            /*chosen_chain=*/cid, /*target_v=*/ch.pending,
                            /*reason=*/"advance:pending_merge");
        ch.pending = -1;
      }
    }
    if (ch.curr + 1 >= static_cast<int>(ch.verts.size())) {
      ch.curr = static_cast<int>(ch.verts.size()) - 2;
    }
    if (ch.curr < 0) ch.curr = 0;
    if (ch.curr != old_curr) {
      ch.cache_stamp = -1;
      // Update current edge parameters for the new curr.
      const int a = ch.verts[ch.curr];
      const int b = ch.verts[ch.curr + 1];
      ch.edge_a = a;
      ch.edge_b = b;
      const auto& p = pts[a];
      const auto& q = pts[b];
      if (std::abs(p.y - q.y) < detail::kEpsMath) {
        ch.edge_slope = 0.0;
        ch.edge_base_x = std::min(p.x, q.x);
      } else {
        const double s = (q.x - p.x) / (p.y - q.y);
        ch.edge_slope = s;
        ch.edge_base_x = p.x + p.y * s;
      }
    }
  }

  double chain_x_at(int cid, const std::vector<Point>& pts) {
    Chain& ch = chains_[cid];
    if (ch.cache_stamp == sweep_stamp_) return ch.cache_x;
    // Only advance when we've swept below the next vertex on this chain.
    if (ch.curr + 1 < static_cast<int>(ch.verts.size()) &&
        pts[ch.verts[ch.curr + 1]].y > sweep_y_eps_) {
      advance_chain(cid, pts);
    }
    const double x = (ch.edge_base_x - sweep_y_ * ch.edge_slope);
    ch.cache_x = x;
    ch.cache_stamp = sweep_stamp_;
    return x;
  }

  int chain_slab_entry_vertex(int cid, const std::vector<Point>& pts) {
    advance_chain(cid, pts);
    return chains_[cid].verts[chains_[cid].curr];
  }

  bool chain_less(int a, int b, const std::vector<Point>& pts) {
    if (a == b) return false;
    const double xa = chain_x_at(a, pts);
    const double xb = chain_x_at(b, pts);
    // IMPORTANT: Must be a strict weak ordering. Do not use epsilon thresholds here.
    if (xa < xb) return true;
    if (xb < xa) return false;
    return a < b;  // stable tie-breaker
  }

  StatusNode* rotate_merge(StatusNode* a, StatusNode* b) {
    if (!a) return b;
    if (!b) return a;
    if (a->prio < b->prio) {
      a->right = rotate_merge(a->right, b);
      return a;
    }
    b->left = rotate_merge(a, b->left);
    return b;
  }

  void split_by_chain(StatusNode* root, int pivot_chain, const std::vector<Point>& pts,
                      StatusNode*& a, StatusNode*& b) {
    // a: chains < pivot_chain, b: chains >= pivot_chain
    if (!root) {
      a = b = nullptr;
      return;
    }
    if (chain_less(root->chain_id, pivot_chain, pts)) {
      split_by_chain(root->right, pivot_chain, pts, root->right, b);
      a = root;
    } else {
      split_by_chain(root->left, pivot_chain, pts, a, root->left);
      b = root;
    }
  }

  void treap_insert(StatusNode*& root, StatusNode* node, const std::vector<Point>& pts) {
    if (!root) {
      root = node;
      return;
    }
    if (node->prio < root->prio) {
      split_by_chain(root, node->chain_id, pts, node->left, node->right);
      root = node;
      return;
    }
    if (chain_less(node->chain_id, root->chain_id, pts)) {
      treap_insert(root->left, node, pts);
    } else {
      treap_insert(root->right, node, pts);
    }
  }

  void treap_erase(StatusNode*& root, int chain_id, const std::vector<Point>& pts) {
    if (!root) return;
    if (root->chain_id == chain_id) {
      StatusNode* old = root;
      root = rotate_merge(root->left, root->right);
      // Nodes are stored in a pool; do not delete.
      (void)old;
      return;
    }
    if (chain_less(chain_id, root->chain_id, pts)) {
      treap_erase(root->left, chain_id, pts);
    } else {
      treap_erase(root->right, chain_id, pts);
    }
  }

  // Predecessor chain strictly left of x_query.
  // Return chain whose x_at(sweep_y_) is maximal strictly less than x_query.
  int predecessor_chain_by_x(StatusNode* root, const std::vector<Point>& pts, double x_query) {
    int best = -1;
    while (root) {
      const double xcur = chain_x_at(root->chain_id, pts);
      if (xcur < x_query) {
        best = root->chain_id;
        root = root->right;
      } else {
        root = root->left;
      }
    }
    return best;
  }

  void dump_status(std::ostream& os, StatusNode* root, const std::vector<Point>& pts) {
    if (!root) return;
    dump_status(os, root->left, pts);
    const int cid = root->chain_id;
    const double x = chain_x_at(cid, pts);
    os << "  chain_id=" << cid << " x=" << x
       << " upper=" << chains_[cid].upper << " lower=" << chains_[cid].lower
       << " pending=" << chains_[cid].pending << "\n";
    dump_status(os, root->right, pts);
  }

  StatusNode* find_node(StatusNode* root, int chain_id, const std::vector<Point>& pts) {
    while (root) {
      if (root->chain_id == chain_id) return root;
      if (chain_less(chain_id, root->chain_id, pts)) root = root->left;
      else root = root->right;
    }
    return nullptr;
  }

  StatusNode* max_node(StatusNode* root) {
    if (!root) return nullptr;
    while (root->right) root = root->right;
    return root;
  }

  StatusNode* min_node(StatusNode* root) {
    if (!root) return nullptr;
    while (root->left) root = root->left;
    return root;
  }

  // NOTE: We intentionally do NOT rely on parent pointers for predecessor/successor.
  // Parent pointers are easy to invalidate in recursive treap operations; relying on them
  // can yield "predecessor not found" even when the node exists.
  int predecessor_of_chain(StatusNode* root, int chain_id, const std::vector<Point>& pts) {
    int best = -1;
    while (root) {
      if (root->chain_id == chain_id) {
        // predecessor is max in left subtree if it exists
        if (root->left) {
          StatusNode* m = max_node(root->left);
          return m ? m->chain_id : -1;
        }
        return best;
      }
      if (chain_less(root->chain_id, chain_id, pts)) {
        best = root->chain_id;
        root = root->right;
      } else {
        root = root->left;
      }
    }
    return -1;
  }

  int successor_of_chain(StatusNode* root, int chain_id, const std::vector<Point>& pts) {
    int best = -1;
    while (root) {
      if (root->chain_id == chain_id) {
        if (root->right) {
          StatusNode* m = min_node(root->right);
          return m ? m->chain_id : -1;
        }
        return best;
      }
      if (chain_less(chain_id, root->chain_id, pts)) {
        best = root->chain_id;
        root = root->left;
      } else {
        root = root->right;
      }
    }
    return -1;
  }

  void free_status(StatusNode* /*root*/) {
    // no-op: status nodes are stored in a pool
  }

  // --- Phase 1: vertex types + chains (left-boundary) ---
  bool is_local_max(const std::vector<Point>& pts, int i) const {
    const int n = static_cast<int>(pts.size());
    const int p = (i - 1 + n) % n;
    const int nx = (i + 1) % n;
    return detail::below(pts, p, i) && detail::below(pts, nx, i);
  }

  bool is_local_min(const std::vector<Point>& pts, int i) const {
    const int n = static_cast<int>(pts.size());
    const int p = (i - 1 + n) % n;
    const int nx = (i + 1) % n;
    return detail::above(pts, p, i) && detail::above(pts, nx, i);
  }

  void build_vertex_types_and_chains(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());
    types_.assign(n, VType::Regular);
    max_to_chain_.assign(n, -1);
    min_to_chain_.assign(n, -1);
    chains_.clear();

    // Vertex type classification (paper definition; assumes CCW).
    for (int i = 0; i < n; ++i) {
      const int p = (i - 1 + n) % n;
      const int nx = (i + 1) % n;
      const bool p_below = detail::below(pts, p, i);
      const bool nx_below = detail::below(pts, nx, i);
      // Reuse convexity flag computed in the reflex-count pass.
      const bool convex = (i < static_cast<int>(is_convex_buf_.size())) && (is_convex_buf_[i] != 0);
      if (p_below && nx_below) {
        types_[i] = convex ? VType::Start : VType::Split;
      } else if (!p_below && !nx_below) {
        types_[i] = convex ? VType::End : VType::Merge;
      } else {
        types_[i] = VType::Regular;
      }
    }

    // Build monotone chains by traversing from each local maximum to the next local minimum
    // along the CCW boundary direction.
    //
    // NOTE: This corresponds to the boundary side where the polygon interior lies to the LEFT
    // when traversing downward (upper -> lower). This matches the edge convention in the
    // classical monotone-decomposition sweep (e_i = (v_i, v_{i+1})).
    for (int max_v = 0; max_v < n; ++max_v) {
      if (!(types_[max_v] == VType::Start || types_[max_v] == VType::Split)) continue;
      Chain ch;
      ch.verts.clear();
      ch.verts.push_back(max_v);

      int cur = (max_v + 1) % n;  // CCW traversal
      int safe = 0;
      while (!(types_[cur] == VType::End || types_[cur] == VType::Merge)) {
        ch.verts.push_back(cur);
        cur = (cur + 1) % n;
        if (++safe > n + 5) throw std::runtime_error("chain construction failed (loop)");
      }
      ch.verts.push_back(cur);

      ch.upper = max_v;
      ch.lower = cur;
      ch.reset_runtime();

      const int id = static_cast<int>(chains_.size());
      chains_.push_back(std::move(ch));
      max_to_chain_[max_v] = id;
      if (min_to_chain_[cur] != -1) {
        // unique left-boundary chain per local minimum
        throw std::runtime_error("duplicate left-boundary chain for a local minimum");
      }
      min_to_chain_[cur] = id;
    }
  }

  // --- Phase 2: chain-based monotone decomposition ---
  struct Event {
    int v;
    double y;
    double x;
    VType type;
  };

  void decompose_into_monotone(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());
    
    // Store pts reference for diagonal validation.
    diag_check_pts_ = &pts;
    
    // Reset diagonal state for this decomposition.
    diagonals_.clear();
    // The monotone partition adds O(#split + #merge) diagonals, and both split/merge
    // are reflex vertices, so O(reflex_count_) is a good capacity estimate.
    diag_set_.reset(static_cast<std::size_t>(reflex_count_) + 32);

    // Build event list: local extrema only (Start/Split/End/Merge), as in the paper.
    events_.clear();
    events_.reserve(n);

    for (int i = 0; i < n; ++i) {
      if (types_[i] == VType::Start || types_[i] == VType::Split ||
          types_[i] == VType::End || types_[i] == VType::Merge) {
        events_.push_back({i, pts[i].y, pts[i].x, types_[i]});
      }
    }
    std::sort(events_.begin(), events_.end(), [](const Event& a, const Event& b) {
      if (a.y > b.y) return true;
      if (a.y < b.y) return false;
      if (a.x < b.x) return true;
      if (a.x > b.x) return false;
      return a.v < b.v;
    });

    // Reset runtime chain state + initialize edge parameters for curr=0.
    for (auto& c : chains_) {
      c.reset_runtime();
      const int a = c.verts[0];
      const int b = c.verts[1];
      c.edge_a = a;
      c.edge_b = b;
      const auto& p = pts[a];
      const auto& q = pts[b];
      if (std::abs(p.y - q.y) < detail::kEpsMath) {
        c.edge_slope = 0.0;
        c.edge_base_x = std::min(p.x, q.x);
      } else {
        const double s = (q.x - p.x) / (p.y - q.y);
        c.edge_slope = s;
        c.edge_base_x = p.x + p.y * s;
      }
    }
    diag_records_.clear();
    pending_records_.clear();

    // Status structure: randomized treap ordered by chain_x_at(sweep_y_) (strict).
    StatusNode* status = nullptr;

    const int chain_count = static_cast<int>(chains_.size());
    if (static_cast<int>(status_pool_.size()) != chain_count) {
      status_pool_.resize(chain_count);
    }

    auto insert_chain = [&](int cid) {
      StatusNode* node = &status_pool_[cid];
      node->chain_id = cid;
      node->prio = next_prio();
      node->left = node->right = nullptr;
      treap_insert(status, node, pts);
    };

    auto remove_chain = [&](int cid) {
      treap_erase(status, cid, pts);
    };

    for (const auto& ev : events_) {
      ++sweep_stamp_;
      const int v = ev.v;
      const double vy = pts[v].y;
      // Sweep convention:
      // - For START/SPLIT (local maxima), we query the status just *below* v so the
      //   outgoing edges from v are considered active.
      // - For END/MERGE (local minima), we query the status just *above* v so the
      //   incoming edges ending at v are still active (and we don't extrapolate
      //   x_at(y) past a chain's lower endpoint).
      if (ev.type == VType::End || ev.type == VType::Merge) {
        sweep_y_ = detail::nextafter_up(vy);
      } else {
        sweep_y_ = detail::nextafter_down(vy);
      }
      sweep_y_eps_ = sweep_y_ + detail::kEpsGeom;

      if (ev.type == VType::Start) {
        const int cid = max_to_chain_[v];
        if (cid < 0) throw std::runtime_error("missing left chain for Start");
        insert_chain(cid);
        chains_[cid].pending = -1;
      } else if (ev.type == VType::End) {
        const int rid = min_to_chain_[v];
        if (rid < 0) throw std::runtime_error("missing left chain for End");
        advance_chain(rid, pts);
        if (chains_[rid].pending != -1) {
          add_diagonal_unique(v, chains_[rid].pending, /*event_v=*/v, ev.type,
                              /*chosen_chain=*/rid, /*target_v=*/chains_[rid].pending,
                              /*reason=*/"end:pending_merge");
          chains_[rid].pending = -1;
        }
        remove_chain(rid);
      } else if (ev.type == VType::Split) {
        const int left_id = predecessor_chain_by_x(status, pts, pts[v].x);
        if (left_id >= 0) {
          int tgt = -1;
          const char* reason = "split:slab_entry";
          if (chains_[left_id].pending != -1) {
            tgt = chains_[left_id].pending;
            chains_[left_id].pending = -1;
            reason = "split:pending_merge";
          } else {
            tgt = chain_slab_entry_vertex(left_id, pts);
          }
          add_diagonal_unique(v, tgt, /*event_v=*/v, ev.type,
                              /*chosen_chain=*/left_id, /*target_v=*/tgt,
                              /*reason=*/reason);
        }
        const int cid = max_to_chain_[v];
        if (cid < 0) throw std::runtime_error("missing left chain for Split");
        insert_chain(cid);
        chains_[cid].pending = -1;
      } else if (ev.type == VType::Merge) {
        const int rid = min_to_chain_[v];
        if (rid < 0) throw std::runtime_error("missing left chain for Merge");
        advance_chain(rid, pts);

        if (chains_[rid].pending != -1) {
          add_diagonal_unique(v, chains_[rid].pending, /*event_v=*/v, ev.type,
                              /*chosen_chain=*/rid, /*target_v=*/chains_[rid].pending,
                              /*reason=*/"merge:term_pending_merge");
          chains_[rid].pending = -1;
        }

        remove_chain(rid);

        // After removing R, find chain immediately left of v.
        const int left_id = predecessor_chain_by_x(status, pts, pts[v].x);
        if (left_id >= 0) {
          if (chains_[left_id].pending != -1) {
            add_diagonal_unique(v, chains_[left_id].pending, /*event_v=*/v, ev.type,
                                /*chosen_chain=*/left_id, /*target_v=*/chains_[left_id].pending,
                                /*reason=*/"merge:left_pending_merge");
            chains_[left_id].pending = -1;
          }
          chains_[left_id].pending = v;
#ifndef NDEBUG
          pending_records_.push_back({left_id, v, v, "merge:set_pending"});
#endif
        }
      } else {
        // no other event types
      }
    }

    free_status(status);
  }

  // --- Phase 3: face extraction + monotone triangulation ---
  void triangulate_monotone_faces(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());

    // Build adjacency in a flat CSR-like layout to avoid O(n) small allocations
    // from vector<vector<int>>. This significantly reduces Phase 3 overhead.
    adj_deg_.assign(n, 2);
    for (const auto& [a, b] : diagonals_) {
      ++adj_deg_[a];
      ++adj_deg_[b];
    }

    adj_off_.assign(n + 1, 0);
    for (int u = 0; u < n; ++u) adj_off_[u + 1] = adj_off_[u] + adj_deg_[u];
    const int m_dir = adj_off_[n];
    dbg_half_edges_ = m_dir;

    adj_flat_.assign(m_dir, -1);
    adj_cur_.resize(n + 1);
    std::copy(adj_off_.begin(), adj_off_.end(), adj_cur_.begin());

    // Polygon boundary edges.
    for (int i = 0; i < n; ++i) {
      const int j = (i + 1) % n;
      adj_flat_[adj_cur_[i]++] = j;
      adj_flat_[adj_cur_[j]++] = i;
    }
    // Decomposition diagonals.
    for (const auto& [a, b] : diagonals_) {
      adj_flat_[adj_cur_[a]++] = b;
      adj_flat_[adj_cur_[b]++] = a;
    }

    // Order neighbors around each vertex in the outerplanar embedding induced by the
    // polygon boundary order (0..n-1). For a non-crossing diagonal set, this provides
    // the correct planar rotation system without any geometric angle computation.
    //
    // For performance: per-vertex sort for smaller n (lower constant), and a linear-time
    // radix counting sort for larger n.
    static constexpr int kSmallPhase3N = 20000;
    const bool small_phase3 = (n <= kSmallPhase3N);
    if (small_phase3) {
      auto key = [&](int u, int v) -> int {
        int d = v - u;
        if (d <= 0) d += n;
        return d;  // 1..n-1, increasing along CCW boundary order
      };
      for (int u = 0; u < n; ++u) {
        const int beg = adj_off_[u];
        const int end = adj_off_[u + 1];
        if (end - beg <= 2) continue;
        std::sort(adj_flat_.begin() + beg, adj_flat_.begin() + end,
                  [&](int a, int b) { return key(u, a) < key(u, b); });
      }
    } else {
      // Linear-time neighbor ordering via radix counting sort over directed half-edges.
      if (static_cast<int>(halfedge_src_.size()) != m_dir) halfedge_src_.resize(m_dir);
      for (int u = 0; u < n; ++u) {
        for (int e = adj_off_[u]; e < adj_off_[u + 1]; ++e) halfedge_src_[e] = u;
      }
      if (static_cast<int>(halfedge_key_.size()) != m_dir) halfedge_key_.resize(m_dir);
      for (int e = 0; e < m_dir; ++e) {
        const int u = halfedge_src_[e];
        const int v = adj_flat_[e];
        int d = v - u;
        if (d <= 0) d += n;
        halfedge_key_[e] = d;
      }

      halfedge_idx1_.resize(m_dir);
      halfedge_idx2_.resize(m_dir);
      for (int i = 0; i < m_dir; ++i) halfedge_idx1_[i] = i;

      // Pass 1: stable sort by key.
      count_buf_.assign(n + 1, 0);
      for (int i = 0; i < m_dir; ++i) ++count_buf_[halfedge_key_[i]];
      int sum = 0;
      for (int k = 0; k <= n; ++k) {
        const int c = count_buf_[k];
        count_buf_[k] = sum;
        sum += c;
      }
      for (int i = 0; i < m_dir; ++i) {
        const int e = halfedge_idx1_[i];
        const int k = halfedge_key_[e];
        halfedge_idx2_[count_buf_[k]++] = e;
      }

      // Pass 2: stable sort by src.
      count_buf_.assign(n + 1, 0);
      for (int i = 0; i < m_dir; ++i) ++count_buf_[halfedge_src_[halfedge_idx2_[i]]];
      sum = 0;
      for (int u = 0; u <= n; ++u) {
        const int c = count_buf_[u];
        count_buf_[u] = sum;
        sum += c;
      }
      for (int i = 0; i < m_dir; ++i) {
        const int e = halfedge_idx2_[i];
        const int u = halfedge_src_[e];
        halfedge_idx1_[count_buf_[u]++] = e;
      }

      adj_flat_tmp_.resize(m_dir);
      for (int i = 0; i < m_dir; ++i) adj_flat_tmp_[i] = adj_flat_[halfedge_idx1_[i]];
      adj_flat_.swap(adj_flat_tmp_);

      // Rebuild src array for the reordered adjacency.
      if (static_cast<int>(halfedge_src_.size()) != m_dir) halfedge_src_.resize(m_dir);
      for (int u = 0; u < n; ++u) {
        for (int e = adj_off_[u]; e < adj_off_[u + 1]; ++e) halfedge_src_[e] = u;
      }
    }

    // Face extraction: traverse directed half-edges in the planar embedding.
    // Reverse half-edge map: for each directed half-edge e=(u->v), store rev(e)=(v->u)
    // as a global index in adj_flat_. Computed in O(n) time via radix counting sort on
    // the undirected key (min(u,v), max(u,v)).
    if (static_cast<int>(halfedge_src_.size()) != m_dir) halfedge_src_.resize(m_dir);
    for (int u = 0; u < n; ++u) {
      for (int e = adj_off_[u]; e < adj_off_[u + 1]; ++e) halfedge_src_[e] = u;
    }

    halfedge_rev_e_.assign(m_dir, -1);
    halfedge_idx1_.resize(m_dir);
    halfedge_idx2_.resize(m_dir);
    for (int i = 0; i < m_dir; ++i) halfedge_idx1_[i] = i;

    // Pass 1: stable sort by b = max(u, v).
    count_buf_.assign(n + 1, 0);
    for (int e = 0; e < m_dir; ++e) {
      const int u = halfedge_src_[e];
      const int v = adj_flat_[e];
      const int b = (u > v) ? u : v;
      ++count_buf_[b];
    }
    int sum2 = 0;
    for (int b = 0; b <= n; ++b) {
      const int c = count_buf_[b];
      count_buf_[b] = sum2;
      sum2 += c;
    }
    for (int i = 0; i < m_dir; ++i) {
      const int e = halfedge_idx1_[i];
      const int u = halfedge_src_[e];
      const int v = adj_flat_[e];
      const int b = (u > v) ? u : v;
      halfedge_idx2_[count_buf_[b]++] = e;
    }

    // Pass 2: stable sort by a = min(u, v).
    count_buf_.assign(n + 1, 0);
    for (int i = 0; i < m_dir; ++i) {
      const int e = halfedge_idx2_[i];
      const int u = halfedge_src_[e];
      const int v = adj_flat_[e];
      const int a = (u < v) ? u : v;
      ++count_buf_[a];
    }
    sum2 = 0;
    for (int a = 0; a <= n; ++a) {
      const int c = count_buf_[a];
      count_buf_[a] = sum2;
      sum2 += c;
    }
    for (int i = 0; i < m_dir; ++i) {
      const int e = halfedge_idx2_[i];
      const int u = halfedge_src_[e];
      const int v = adj_flat_[e];
      const int a = (u < v) ? u : v;
      halfedge_idx1_[count_buf_[a]++] = e;
    }

    for (int i = 0; i + 1 < m_dir; i += 2) {
      const int e1 = halfedge_idx1_[i];
      const int e2 = halfedge_idx1_[i + 1];
      halfedge_rev_e_[e1] = e2;
      halfedge_rev_e_[e2] = e1;
    }

    // Validate reverse half-edge pairing (debug-only; can be costly at scale).
#ifndef NDEBUG
    {
      int bad_pairs = 0;
      for (int e = 0; e < m_dir; ++e) {
        const int rev_e = halfedge_rev_e_[e];
        if (rev_e < 0 || rev_e >= m_dir) { ++bad_pairs; continue; }
        const int src_e = halfedge_src_[e];
        const int dst_e = adj_flat_[e];
        const int src_rev = halfedge_src_[rev_e];
        const int dst_rev = adj_flat_[rev_e];
        if (src_e != dst_rev || dst_e != src_rev) ++bad_pairs;
      }
      if (bad_pairs > 0) {
        throw std::runtime_error("bad half-edge pairing: " + std::to_string(bad_pairs));
      }
    }
#endif

    // Mark the outer face half-edges so we don't triangulate it.
    halfedge_used_.assign(m_dir, 0);
    int used_count = 0;
    int outer_start_e = -1;
    if (n >= 2) {
      const int u0 = 0;
      const int v0 = n - 1;
      const int beg0 = adj_off_[u0];
      const int end0 = adj_off_[u0 + 1];
      if (end0 > beg0 && adj_flat_[end0 - 1] == v0) {
        outer_start_e = end0 - 1;
      } else {
        for (int e = beg0; e < end0; ++e) {
          if (adj_flat_[e] == v0) { outer_start_e = e; break; }
        }
      }
    }
    if (outer_start_e >= 0) {
      const int start_u = 0;
      int prev_v = start_u;
      int curr = adj_flat_[outer_start_e];
      halfedge_used_[outer_start_e] = 1;
      ++used_count;
      int idx_prev = halfedge_rev_e_[outer_start_e] - adj_off_[curr];
      int iterations = 0;
      while (curr != start_u && iterations < 2 * m_dir) {
        const int deg_curr = adj_off_[curr + 1] - adj_off_[curr];
        if (idx_prev < 0 || idx_prev >= deg_curr) break;
        const int next_i = (idx_prev - 1 + deg_curr) % deg_curr;
        const int e2 = adj_off_[curr] + next_i;
        const int next_v = adj_flat_[e2];
        if (!halfedge_used_[e2]) {
          halfedge_used_[e2] = 1;
          ++used_count;
        }
        const int rev_e2 = halfedge_rev_e_[e2];
        prev_v = curr;
        curr = next_v;
        idx_prev = rev_e2 - adj_off_[next_v];
        ++iterations;
      }
    }

#ifndef NDEBUG
    dbg_faces_ = 0;
    dbg_failed_walks_ = 0;
    dbg_degenerate_faces_ = 0;
    dbg_iter_limit_faces_ = 0;
    dbg_cw_faces_ = 0;
    dbg_sum_face_verts_ = 0;
#endif

    face_buf_.clear();
    face_buf_.reserve(32);

    for (int start_u = 0; start_u < n; ++start_u) {
      const int deg_u = adj_off_[start_u + 1] - adj_off_[start_u];
      for (int start_i = 0; start_i < deg_u; ++start_i) {
        const int start_e = adj_off_[start_u] + start_i;
        if (halfedge_used_[start_e]) continue;

        face_buf_.clear();
        face_buf_.push_back(start_u);

        int prev_v = start_u;
        int curr = adj_flat_[start_e];

        halfedge_used_[start_e] = 1;
        ++used_count;
        int idx_prev = halfedge_rev_e_[start_e] - adj_off_[curr];

        double area2 = 0.0;
        int iterations = 0;
        while (curr != start_u && iterations < 2 * m_dir) {
          face_buf_.push_back(curr);
          area2 += pts[prev_v].x * pts[curr].y - pts[curr].x * pts[prev_v].y;

          const int deg_curr = adj_off_[curr + 1] - adj_off_[curr];
          if (idx_prev < 0 || idx_prev >= deg_curr) break;

          const int next_i = (idx_prev - 1 + deg_curr) % deg_curr;
          const int e2 = adj_off_[curr] + next_i;
          const int next_v = adj_flat_[e2];

          if (!halfedge_used_[e2]) {
            halfedge_used_[e2] = 1;
            ++used_count;
          }

          const int rev_e2 = halfedge_rev_e_[e2];
          prev_v = curr;
          curr = next_v;
          idx_prev = rev_e2 - adj_off_[next_v];
          ++iterations;
        }

        if (curr == start_u && face_buf_.size() >= 3 && iterations < 2 * m_dir) {
          area2 += pts[prev_v].x * pts[start_u].y - pts[start_u].x * pts[prev_v].y;
          // Ensure CCW orientation for monotone triangulation.
          if (area2 < -1e-12) std::reverse(face_buf_.begin(), face_buf_.end());
#ifndef NDEBUG
          ++dbg_faces_;
          dbg_sum_face_verts_ += static_cast<long long>(face_buf_.size());
#endif
          triangulate_one_monotone_face(pts, face_buf_);
        }
#ifndef NDEBUG
        else if (curr != start_u) {
          ++dbg_failed_walks_;
        } else if (face_buf_.size() < 3) {
          ++dbg_degenerate_faces_;
        } else if (iterations >= 2 * m_dir) {
          ++dbg_iter_limit_faces_;
        }
#endif
      }
    }

    dbg_used_half_edges_ = used_count;

#ifndef NDEBUG
    // Count unused half-edges + Euler mismatch (debug-only).
    int unused = 0;
    for (int e = 0; e < m_dir; ++e) {
      if (!halfedge_used_[e]) ++unused;
    }
    dbg_unused_halfedges_ = unused;
    dbg_total_walks_ = dbg_faces_ + dbg_cw_faces_ + dbg_failed_walks_ + dbg_degenerate_faces_ + dbg_iter_limit_faces_;
    dbg_euler_mismatch_ = (2 + static_cast<int>(diagonals_.size())) - (dbg_faces_ + dbg_cw_faces_);
#endif
  }

  // Check if segment (a,b) properly intersects segment (c,d).
  // Returns true if they cross in their interiors (excluding shared endpoints).
  static bool segments_cross(const Point& a, const Point& b, const Point& c, const Point& d) {
    const double o1 = detail::cross(a, b, c);
    const double o2 = detail::cross(a, b, d);
    const double o3 = detail::cross(c, d, a);
    const double o4 = detail::cross(c, d, b);
    // Proper intersection: points on opposite sides of each line.
    const double eps = detail::kEpsGeom;
    if ((o1 > eps && o2 < -eps) || (o1 < -eps && o2 > eps)) {
      if ((o3 > eps && o4 < -eps) || (o3 < -eps && o4 > eps)) {
        return true;
      }
    }
    return false;
  }

  // Check if diagonal from face[i] to face[j] is valid for ear clipping.
  // The diagonal is valid if:
  // 1. It doesn't cross any face boundary edge
  // 2. It lies inside the polygon (midpoint is inside)
  bool is_valid_ear_diagonal(const std::vector<Point>& pts, const std::vector<int>& face_list,
                              int i, int j) const {
    const int m = static_cast<int>(face_list.size());
    const int vi = face_list[i];
    const int vj = face_list[j];
    const Point& pi = pts[vi];
    const Point& pj = pts[vj];

    // Check crossing with all non-adjacent edges.
    for (int k = 0; k < m; ++k) {
      const int next_k = (k + 1) % m;
      if (k == i || k == j || next_k == i || next_k == j) continue;
      if (segments_cross(pi, pj, pts[face_list[k]], pts[face_list[next_k]])) {
        return false;
      }
    }
    return true;
  }

  // Triangulate a face using ear clipping (works for any simple polygon).
  // O(m^2) time, but faces are typically small.
  void triangulate_face_ear_clipping(const std::vector<Point>& pts, const std::vector<int>& face) {
    const int m = static_cast<int>(face.size());
    if (m < 3) return;
    if (m == 3) {
      triangles_.push_back({face[0], face[1], face[2]});
      return;
    }

    // Work with a mutable list of remaining vertices.
    ear_clip_list_.clear();
    ear_clip_list_.reserve(m);
    for (int i = 0; i < m; ++i) {
      ear_clip_list_.push_back(face[i]);
    }

    // Compute signed area to determine winding (CCW = positive).
    double area2 = 0.0;
    for (int i = 0; i < m; ++i) {
      const int j = (i + 1) % m;
      area2 += pts[face[i]].x * pts[face[j]].y - pts[face[j]].x * pts[face[i]].y;
    }
    const bool ccw = area2 > 0;

    // Use a relaxed epsilon for ear tests to handle numerical edge cases.
    const double ear_eps = 1e-14;

    int iterations = 0;
    const int max_iterations = m * m * 2;  // Safety limit.

    while (ear_clip_list_.size() > 3 && iterations < max_iterations) {
      const int n = static_cast<int>(ear_clip_list_.size());
      bool found_ear = false;
      int best_ear = -1;
      double best_convexity = 0.0;

      // First pass: find the most convex ear.
      for (int i = 0; i < n; ++i) {
        const int prev_i = (i - 1 + n) % n;
        const int next_i = (i + 1) % n;

        const int vprev = ear_clip_list_[prev_i];
        const int vcurr = ear_clip_list_[i];
        const int vnext = ear_clip_list_[next_i];

        // Check if this is a convex vertex (ear candidate).
        const double c = detail::cross(pts[vprev], pts[vcurr], pts[vnext]);
        const bool is_convex = ccw ? (c > -ear_eps) : (c < ear_eps);
        if (!is_convex) continue;

        // Check if any other vertex is inside the triangle.
        bool has_point_inside = false;
        for (int k = 0; k < n; ++k) {
          if (k == prev_i || k == i || k == next_i) continue;
          const int vk = ear_clip_list_[k];
          // Check if vk is strictly inside triangle (vprev, vcurr, vnext).
          const double c1 = detail::cross(pts[vprev], pts[vcurr], pts[vk]);
          const double c2 = detail::cross(pts[vcurr], pts[vnext], pts[vk]);
          const double c3 = detail::cross(pts[vnext], pts[vprev], pts[vk]);
          bool inside;
          if (ccw) {
            inside = c1 > ear_eps && c2 > ear_eps && c3 > ear_eps;
          } else {
            inside = c1 < -ear_eps && c2 < -ear_eps && c3 < -ear_eps;
          }
          if (inside) {
            has_point_inside = true;
            break;
          }
        }
        if (has_point_inside) continue;

        // This is a valid ear. Track the most convex one.
        const double convexity = ccw ? c : -c;
        if (best_ear < 0 || convexity > best_convexity) {
          best_ear = i;
          best_convexity = convexity;
        }
      }

      if (best_ear >= 0) {
        const int i = best_ear;
        const int prev_i = (i - 1 + n) % n;
        const int next_i = (i + 1) % n;
        triangles_.push_back({ear_clip_list_[prev_i], ear_clip_list_[i], ear_clip_list_[next_i]});
        ear_clip_list_.erase(ear_clip_list_.begin() + i);
        found_ear = true;
      }

      if (!found_ear) {
        // No valid ear found. This can happen with degenerate polygons.
        // Find the vertex closest to collinear and remove it.
        double min_abs_cross = std::numeric_limits<double>::max();
        int min_idx = 0;
        for (int i = 0; i < n; ++i) {
          const int prev_i = (i - 1 + n) % n;
          const int next_i = (i + 1) % n;
          const double c = std::abs(detail::cross(pts[ear_clip_list_[prev_i]],
                                                   pts[ear_clip_list_[i]],
                                                   pts[ear_clip_list_[next_i]]));
          if (c < min_abs_cross) {
            min_abs_cross = c;
            min_idx = i;
          }
        }
        const int i = min_idx;
        const int prev_i = (i - 1 + n) % n;
        const int next_i = (i + 1) % n;
        triangles_.push_back({ear_clip_list_[prev_i], ear_clip_list_[i], ear_clip_list_[next_i]});
        ear_clip_list_.erase(ear_clip_list_.begin() + i);
      }

      ++iterations;
    }

    // Handle remaining triangle.
    if (ear_clip_list_.size() == 3) {
      triangles_.push_back({ear_clip_list_[0], ear_clip_list_[1], ear_clip_list_[2]});
    }
  }

  void triangulate_one_monotone_face(const std::vector<Point>& pts, const std::vector<int>& face) {
    const int m = static_cast<int>(face.size());
    if (m < 3) return;
    if (m == 3) {
      triangles_.push_back({face[0], face[1], face[2]});
      return;
    }

    // Triangulate a y-monotone face in O(m) time (PolyPartition-style, robust).
    // Work in the index space of `face` (0..m-1) and emit triangles in original vertex ids.
    int topindex = 0;
    int bottomindex = 0;
    for (int i = 1; i < m; ++i) {
      if (detail::below(pts, face[i], face[bottomindex])) bottomindex = i;
      if (detail::below(pts, face[topindex], face[i])) topindex = i;
    }

    // Verify monotonicity only in debug builds (performance-critical in benchmarks).
#ifndef NDEBUG
    {
      int i = topindex;
      while (i != bottomindex) {
        int i2 = i + 1;
        if (i2 >= m) i2 = 0;
        if (!detail::below(pts, face[i2], face[i])) {
          triangulate_face_ear_clipping(pts, face);
          return;
        }
        i = i2;
      }
      i = bottomindex;
      while (i != topindex) {
        int i2 = i + 1;
        if (i2 >= m) i2 = 0;
        if (!detail::below(pts, face[i], face[i2])) {
          triangulate_face_ear_clipping(pts, face);
          return;
        }
        i = i2;
      }
    }
#endif

    // Merge left and right vertex chains into `tmp_sorted_` (priority order).
    // Also store vertex types in `tmp_chain1_`:  1=left chain, -1=right chain, 0=top/bottom.
    tmp_sorted_.assign(m, 0);
    tmp_chain1_.assign(m, 0);
    tmp_sorted_[0] = topindex;
    tmp_chain1_[topindex] = 0;
    int leftindex = topindex + 1;
    if (leftindex >= m) leftindex = 0;
    int rightindex = topindex - 1;
    if (rightindex < 0) rightindex = m - 1;

    for (int k = 1; k < m - 1; ++k) {
      if (leftindex == bottomindex) {
        tmp_sorted_[k] = rightindex;
        rightindex--;
        if (rightindex < 0) rightindex = m - 1;
        tmp_chain1_[tmp_sorted_[k]] = -1;
      } else if (rightindex == bottomindex) {
        tmp_sorted_[k] = leftindex;
        leftindex++;
        if (leftindex >= m) leftindex = 0;
        tmp_chain1_[tmp_sorted_[k]] = 1;
      } else {
        // Pick the higher (more topmost) vertex among left/right.
        if (detail::below(pts, face[leftindex], face[rightindex])) {
          tmp_sorted_[k] = rightindex;
          rightindex--;
          if (rightindex < 0) rightindex = m - 1;
          tmp_chain1_[tmp_sorted_[k]] = -1;
        } else {
          tmp_sorted_[k] = leftindex;
          leftindex++;
          if (leftindex >= m) leftindex = 0;
          tmp_chain1_[tmp_sorted_[k]] = 1;
        }
      }
    }
    tmp_sorted_[m - 1] = bottomindex;
    tmp_chain1_[bottomindex] = 0;

    // Stack-based triangulation.
    tmp_stack_.resize(m);
    int stackptr = 0;
    tmp_stack_[0] = tmp_sorted_[0];
    tmp_stack_[1] = tmp_sorted_[1];
    stackptr = 2;

    for (int k = 2; k < m - 1; ++k) {
      const int vpos = tmp_sorted_[k];
      if (tmp_chain1_[vpos] != tmp_chain1_[tmp_stack_[stackptr - 1]]) {
        for (int j = 0; j < stackptr - 1; ++j) {
          const int a = tmp_stack_[j];
          const int b = tmp_stack_[j + 1];
          if (tmp_chain1_[vpos] == 1) {
            triangles_.push_back({face[b], face[a], face[vpos]});
          } else {
            triangles_.push_back({face[a], face[b], face[vpos]});
          }
        }
        tmp_stack_[0] = tmp_sorted_[k - 1];
        tmp_stack_[1] = tmp_sorted_[k];
        stackptr = 2;
      } else {
        stackptr--;
        while (stackptr > 0) {
          const int p_v = face[vpos];
          if (tmp_chain1_[vpos] == 1) {
            // IsConvex(v, stack[stackptr-1], stack[stackptr])
            const int p_a = face[tmp_stack_[stackptr - 1]];
            const int p_b = face[tmp_stack_[stackptr]];
            if (detail::cross(pts[p_v], pts[p_a], pts[p_b]) > detail::kEpsGeom) {
              triangles_.push_back({p_v, p_a, p_b});
              stackptr--;
            } else {
              break;
            }
          } else {
            // IsConvex(v, stack[stackptr], stack[stackptr-1])
            const int p_a = face[tmp_stack_[stackptr]];
            const int p_b = face[tmp_stack_[stackptr - 1]];
            if (detail::cross(pts[p_v], pts[p_a], pts[p_b]) > detail::kEpsGeom) {
              triangles_.push_back({p_v, p_a, p_b});
              stackptr--;
            } else {
              break;
            }
          }
        }
        stackptr++;
        tmp_stack_[stackptr] = vpos;
        stackptr++;
      }
    }

    const int vpos = tmp_sorted_[m - 1];  // bottom
    for (int j = 0; j < stackptr - 1; ++j) {
      const int a = tmp_stack_[j];
      const int b = tmp_stack_[j + 1];
      if (tmp_chain1_[b] == 1) {
        triangles_.push_back({face[a], face[b], face[vpos]});
      } else {
        triangles_.push_back({face[b], face[a], face[vpos]});
      }
    }
  }
};

}  // namespace reflex_tri


