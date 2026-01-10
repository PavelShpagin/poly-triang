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

static constexpr double kEps = 1e-12;

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
    if (c > kEps) {
      if (sign < 0) return false;
      sign = 1;
    } else if (c < -kEps) {
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
  if (pts[a].y < pts[b].y) return true;
  if (pts[a].y > pts[b].y) return false;
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
  return cross(pts[p], pts[i], pts[nx]) < -kEps;
}

}  // namespace detail

class Triangulator {
public:
  const std::vector<Triangle>& triangulate(std::vector<Point>& pts) {
    triangles_.clear();
    diagonals_.clear();
    reflex_count_ = 0;
    dbg_faces_ = 0;
    dbg_sum_face_verts_ = 0;
    dbg_half_edges_ = 0;
    dbg_used_half_edges_ = 0;

    const int n = static_cast<int>(pts.size());
    if (n < 3) return triangles_;

    // Fast path for convex polygons (including CW order): fan triangulation.
    // We gate this check to small n because, on non-convex inputs, a full convexity scan
    // adds an extra O(n) pass.
    if (n <= 2048 && detail::is_convex_polygon(pts)) {
      triangles_.resize(n - 2);
      for (int i = 0; i < n - 2; ++i) triangles_[i] = {0, i + 1, i + 2};
      return triangles_;
    }

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
      if (c > detail::kEps) is_convex_buf_[i] = 1;
      if (c < -detail::kEps) ++reflex_count_;
    }

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
    // Helper of the current active edge (PolyPartition helper(e)).
    // We only need to special-case MERGE helpers for diagonal insertion, but
    // we still must update the helper on right-chain REGULAR and SPLIT vertices
    // to avoid stale MERGE helpers.
    int helper = -1;         // helper vertex index (any type), or -1 if unset
    int pending = -1;        // cached: helper if types_[helper]==Merge, else -1
    int upper = -1;
    int lower = -1;
    double cache_y = std::numeric_limits<double>::quiet_NaN();
    double cache_x = 0.0;
    // Current edge parameters for fast x_at(y):
    // If edge is non-horizontal, x(y) = base_x - slope * y (derived from the edge endpoints).
    int edge_a = -1;
    int edge_b = -1;
    double edge_base_x = 0.0;
    double edge_slope = 0.0;
    bool edge_horizontal = false;

    void reset_runtime() {
      curr = 0;
      helper = -1;
      pending = -1;
      cache_y = std::numeric_limits<double>::quiet_NaN();
      cache_x = 0.0;
      edge_a = -1;
      edge_b = -1;
      edge_base_x = 0.0;
      edge_slope = 0.0;
      edge_horizontal = false;
    }

    void advance(double y, const std::vector<Point>& pts) {
      // Move down the chain until current edge spans y.
      const int old_curr = curr;
      while (curr + 1 < static_cast<int>(verts.size()) &&
             pts[verts[curr + 1]].y > y + detail::kEps) {
        ++curr;
      }
      if (curr + 1 >= static_cast<int>(verts.size())) {
        curr = static_cast<int>(verts.size()) - 2;
      }
      if (curr < 0) curr = 0;
      if (curr != old_curr) {
        cache_y = std::numeric_limits<double>::quiet_NaN();
        edge_a = -1;
      }
    }

    int slab_entry_vertex() const {
      // Upper vertex of current edge.
      return verts[curr];
    }

    // NOTE: We do NOT expose a generic x_at(y) here because, in the full
    // algorithm, advancing a chain must also resolve pending MERGE helpers at
    // skipped REGULAR vertices. That logic lives in the owning Triangulator.
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
  std::vector<std::pair<int, int>> diagonals_;
  // For uniqueness without an O(|D| log |D|) sort pass.
  DiagSet diag_set_;
  std::vector<DiagRecord> diag_records_;
  std::vector<PendingRecord> pending_records_;
  std::vector<Triangle> triangles_;
  int reflex_count_ = 0;
  std::vector<unsigned char> is_convex_buf_;
#ifdef REFLEX_TRI_TIMING
  double phase1_ms_ = 0.0;
  double phase2_ms_ = 0.0;
  double phase3_ms_ = 0.0;
#endif
  // Debug stats from face extraction (helpful to diagnose degeneracy / bugs).
  int dbg_faces_ = 0;
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

  bool add_diagonal_unique(int u, int v, int event_v, VType event_type, int chosen_chain,
                           int target_v, const char* reason) {
    const int n = static_cast<int>(types_.size());
    // Never insert a boundary edge as a "diagonal" (can create duplicate edges in Phase 3).
    if (n > 0) {
      if ((u + 1) % n == v || (v + 1) % n == u) return false;
    }
    if (u == v) return false;
    if (u > v) std::swap(u, v);
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
      // Only the "interior to right" regular case corresponds to helper(e_{i-1})
      // on this chain (PolyPartition's REGULAR/right case):
      // Below(v_i, v_{i-1}).
      const int prev = (vreg - 1 + n) % n;
      if (!detail::below(pts, vreg, prev)) continue;

      // If helper(e_{i-1}) is MERGE, add diagonal to it.
      if (ch.pending != -1) {
        add_diagonal_unique(vreg, ch.pending, /*event_v=*/vreg, VType::Regular,
                            /*chosen_chain=*/cid, /*target_v=*/ch.pending,
                            /*reason=*/"regular:merge_helper");
      }

      // REGULAR/right: helper(e_i) becomes vreg (non-merge).
      ch.helper = vreg;
      ch.pending = -1;
    }
    if (ch.curr + 1 >= static_cast<int>(ch.verts.size())) {
      ch.curr = static_cast<int>(ch.verts.size()) - 2;
    }
    if (ch.curr < 0) ch.curr = 0;
    if (ch.curr != old_curr) {
      ch.cache_y = std::numeric_limits<double>::quiet_NaN();
      ch.edge_a = -1;
    }
  }

  double chain_x_at(int cid, const std::vector<Point>& pts) {
    Chain& ch = chains_[cid];
    if (ch.cache_y == sweep_y_) return ch.cache_x;
    // Only advance when we've swept below the next vertex on this chain.
    if (ch.curr + 1 < static_cast<int>(ch.verts.size()) &&
        pts[ch.verts[ch.curr + 1]].y > sweep_y_eps_) {
      advance_chain(cid, pts);
    }
    const int a = ch.verts[ch.curr];
    const int b = ch.verts[ch.curr + 1];
    if (a != ch.edge_a || b != ch.edge_b) {
      ch.edge_a = a;
      ch.edge_b = b;
      const auto& p = pts[a];
      const auto& q = pts[b];
      if (std::abs(p.y - q.y) < detail::kEps) {
        ch.edge_horizontal = true;
        ch.edge_base_x = std::min(p.x, q.x);
        ch.edge_slope = 0.0;
      } else {
        ch.edge_horizontal = false;
        // x(y) = p.x + (p.y - y) * (q.x - p.x) / (p.y - q.y)
        //      = (p.x + p.y*s) - y*s  where s=(q.x-p.x)/(p.y-q.y)
        const double s = (q.x - p.x) / (p.y - q.y);
        ch.edge_slope = s;
        ch.edge_base_x = p.x + p.y * s;
      }
    }
    const double x = ch.edge_horizontal ? ch.edge_base_x : (ch.edge_base_x - sweep_y_ * ch.edge_slope);
    ch.cache_y = sweep_y_;
    ch.cache_x = x;
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
    // Reset diagonal state for this decomposition.
    diagonals_.clear();
    // The monotone partition adds O(#split + #merge) diagonals, and both split/merge
    // are reflex vertices, so O(reflex_count_) is a good capacity estimate.
    diag_set_.reset(static_cast<std::size_t>(reflex_count_) + 32);

    // Build event list:
    // - All extrema (START/SPLIT/END/MERGE)
    // - PLUS regular vertices on the *right* chain (PolyPartition REGULAR/left case),
    //   because these vertices update helper(e_j) and can clear a pending MERGE helper.
    //
    // The original extrema-only variant can leave MERGE helpers stale and produce
    // invalid diagonals on some polygons (e.g., MERGE-to-MERGE chords crossing the boundary).
    std::vector<Event> events;
    events.reserve(n);

    for (int i = 0; i < n; ++i) {
      if (types_[i] == VType::Start || types_[i] == VType::Split ||
          types_[i] == VType::End || types_[i] == VType::Merge) {
        events.push_back({i, pts[i].y, pts[i].x, types_[i]});
      } else if (types_[i] == VType::Regular) {
        // Regular vertex on the right chain: prev is below v (interior to left).
        const int p = (i - 1 + n) % n;
        if (!detail::below(pts, i, p)) {
          events.push_back({i, pts[i].y, pts[i].x, VType::Regular});
        }
      }
    }
    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
      if (std::abs(a.y - b.y) > detail::kEps) return a.y > b.y;
      return a.x < b.x;
    });

    // Reset runtime chain state.
    for (auto& c : chains_) c.reset_runtime();
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

    for (const auto& ev : events) {
      const int v = ev.v;
      const double vy = pts[v].y;
      // Sweep convention:
      // - For START/SPLIT (local maxima), we query the status just *below* v so the
      //   outgoing edges from v are considered active.
      // - For END/MERGE (local minima), we query the status just *above* v so the
      //   incoming edges ending at v are still active (and we don't extrapolate
      //   x_at(y) past a chain's lower endpoint).
      if (ev.type == VType::End || ev.type == VType::Merge || ev.type == VType::Regular) {
        sweep_y_ = detail::nextafter_up(vy);
      } else {
        sweep_y_ = detail::nextafter_down(vy);
      }
      sweep_y_eps_ = sweep_y_ + detail::kEps;

      if (ev.type == VType::Start) {
        const int cid = max_to_chain_[v];
        if (cid < 0) throw std::runtime_error("missing left chain for Start");
        insert_chain(cid);
        chains_[cid].helper = v;
        chains_[cid].pending = -1;
      } else if (ev.type == VType::End) {
        const int rid = min_to_chain_[v];
        if (rid < 0) throw std::runtime_error("missing left chain for End");
        // Ensure all skipped REGULAR vertices on R are processed before we act on pending.
        advance_chain(rid, pts);
        if (chains_[rid].pending != -1) {
          add_diagonal_unique(v, chains_[rid].pending, /*event_v=*/v, ev.type,
                              /*chosen_chain=*/rid, /*target_v=*/chains_[rid].pending,
                              /*reason=*/"end:merge_helper");
          chains_[rid].pending = -1;
        }
        remove_chain(rid);
      } else if (ev.type == VType::Split) {
        const int left_id = predecessor_chain_by_x(status, pts, pts[v].x);
        if (left_id >= 0) {
          const int tgt = (chains_[left_id].helper != -1)
                              ? chains_[left_id].helper
                              : chain_slab_entry_vertex(left_id, pts);
          add_diagonal_unique(v, tgt, /*event_v=*/v, ev.type,
                              /*chosen_chain=*/left_id, /*target_v=*/tgt,
                              /*reason=*/(chains_[left_id].pending != -1) ? "split:merge_helper" : "split:helper");
          // helper(e_j) becomes v (SPLIT is non-merge)
          chains_[left_id].helper = v;
          chains_[left_id].pending = -1;
        }
        const int cid = max_to_chain_[v];
        if (cid < 0) throw std::runtime_error("missing left chain for Split");
        insert_chain(cid);
        chains_[cid].helper = v;
        chains_[cid].pending = -1;
      } else if (ev.type == VType::Merge) {
        const int rid = min_to_chain_[v];
        if (rid < 0) throw std::runtime_error("missing left chain for Merge");

        // Ensure all skipped REGULAR vertices on R are processed before we act on pending.
        advance_chain(rid, pts);
        if (chains_[rid].pending != -1) {
          add_diagonal_unique(v, chains_[rid].pending, /*event_v=*/v, ev.type,
                              /*chosen_chain=*/rid, /*target_v=*/chains_[rid].pending,
                              /*reason=*/"merge:term_merge_helper");
          chains_[rid].pending = -1;
        }

        // Remove R from the status structure (edge e_{i-1}).
        remove_chain(rid);

        // After removal, find the chain immediately left of v (edge e_j).
        const int left_id = predecessor_chain_by_x(status, pts, pts[v].x);
        if (left_id >= 0) {
          if (chains_[left_id].pending != -1) {
            add_diagonal_unique(v, chains_[left_id].pending, /*event_v=*/v, ev.type,
                                /*chosen_chain=*/left_id, /*target_v=*/chains_[left_id].pending,
                                /*reason=*/"merge:left_merge_helper");
            chains_[left_id].pending = -1;
          }
          // helper(e_j) becomes v (MERGE is a merge helper)
          chains_[left_id].helper = v;
          chains_[left_id].pending = v;
          #ifndef NDEBUG
          pending_records_.push_back({left_id, v, v, "merge:set_pending"});
          #endif
        }
      } else {
        // Regular vertex on the right chain: update helper of left edge e_j.
        // If helper(e_j) is MERGE, add the diagonal; then helper(e_j) becomes v.
        const int left_id = predecessor_chain_by_x(status, pts, pts[v].x);

        if (left_id >= 0) {
          if (chains_[left_id].pending != -1) {
            add_diagonal_unique(v, chains_[left_id].pending, /*event_v=*/v, ev.type,
                                /*chosen_chain=*/left_id, /*target_v=*/chains_[left_id].pending,
                                /*reason=*/"regular:right_merge_helper");
          }
          chains_[left_id].helper = v;
          chains_[left_id].pending = -1;
        }
      }
    }

    free_status(status);
  }

  // --- Phase 3: face extraction + monotone triangulation ---
  void triangulate_monotone_faces(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());

    // Build adjacency in a flat CSR-like layout to avoid O(n) small allocations
    // from vector<vector<int>>. This significantly reduces Phase 3 overhead.
    std::vector<int> deg(n, 2);
    for (const auto& [a, b] : diagonals_) {
      ++deg[a];
      ++deg[b];
    }

    std::vector<int> off(n + 1, 0);
    for (int u = 0; u < n; ++u) off[u + 1] = off[u] + deg[u];
    const int m_dir = off[n];
    dbg_half_edges_ = m_dir;

    std::vector<int> adj(m_dir, -1);
    std::vector<int> cur = off;

    // Polygon boundary edges.
    for (int i = 0; i < n; ++i) {
      const int j = (i + 1) % n;
      adj[cur[i]++] = j;
      adj[cur[j]++] = i;
    }
    // Decomposition diagonals.
    for (const auto& [a, b] : diagonals_) {
      adj[cur[a]++] = b;
      adj[cur[b]++] = a;
    }

    // Sort neighbors around each vertex by angle (CCW), without atan2.
    for (int u = 0; u < n; ++u) {
      const int beg = off[u];
      const int end = off[u + 1];
      if (end - beg <= 2) continue;  // ordering is irrelevant for degree <= 2
      auto half = [&](int v) -> int {
        const double dx = pts[v].x - pts[u].x;
        const double dy = pts[v].y - pts[u].y;
        if (dy > 0) return 0;
        if (dy < 0) return 1;
        return (dx >= 0) ? 0 : 1;
      };
      std::sort(adj.begin() + beg, adj.begin() + end, [&](int a, int b) {
        const int ha = half(a);
        const int hb = half(b);
        if (ha != hb) return ha < hb;
        const double ax = pts[a].x - pts[u].x;
        const double ay = pts[a].y - pts[u].y;
        const double bx = pts[b].x - pts[u].x;
        const double by = pts[b].y - pts[u].y;
        const double cr = ax * by - ay * bx;
        // IMPORTANT: strict weak ordering (no eps threshold here).
        if (cr > 0) return true;
        if (cr < 0) return false;
        const double da = ax * ax + ay * ay;
        const double db = bx * bx + by * by;
        if (da < db) return true;
        if (db < da) return false;
        return a < b;
      });
    }

    // Face extraction: traverse directed half-edges in the planar embedding.
    std::vector<unsigned char> used(m_dir, 0);
    int used_count = 0;

    // Reverse lookup: for directed edge (u -> adj[u][i]), store the index of u in adj[v].
    // In practice, neighbor degrees are small; a linear scan is faster than a global hash map.
    std::vector<int> rev_pos(m_dir, -1);
    for (int u = 0; u < n; ++u) {
      const int beg_u = off[u];
      const int end_u = off[u + 1];
      for (int e = beg_u; e < end_u; ++e) {
        const int v = adj[e];
        const int beg_v = off[v];
        const int end_v = off[v + 1];
        int j = -1;
        for (int t = beg_v; t < end_v; ++t) {
          if (adj[t] == u) { j = t - beg_v; break; }
        }
        rev_pos[e] = j;
      }
    }

    // Triangulate each interior face on-the-fly (no storing of all faces).
    // With CCW neighbor ordering and the "prev neighbor" face-walk rule, interior
    // faces are traversed CCW (positive signed area) while the outer face is CW.
    if (static_cast<int>(is_left_buf_.size()) != n) {
      is_left_buf_.assign(n, 0);
    } else {
      std::fill(is_left_buf_.begin(), is_left_buf_.end(), 0);
    }
    is_left_touched_.clear();
    dbg_faces_ = 0;
    dbg_sum_face_verts_ = 0;

    std::vector<int> face;
    face.reserve(32);

    for (int start_u = 0; start_u < n; ++start_u) {
      const int deg_u = off[start_u + 1] - off[start_u];
      for (int start_i = 0; start_i < deg_u; ++start_i) {
        const int start_e = off[start_u] + start_i;
        if (used[start_e]) continue;

        face.clear();
        face.push_back(start_u);

        int curr = adj[start_e];

        used[start_e] = 1;
        ++used_count;
        int idx_prev = rev_pos[start_e];

        double area2 = 0.0;
        int prev_v = start_u;

        int iterations = 0;
        while (curr != start_u && iterations < 2 * m_dir) {
          face.push_back(curr);
          area2 += pts[prev_v].x * pts[curr].y - pts[curr].x * pts[prev_v].y;

          const int deg_curr = off[curr + 1] - off[curr];
          if (idx_prev < 0 || idx_prev >= deg_curr) break;

          const int next_i = (idx_prev - 1 + deg_curr) % deg_curr;
          const int e2 = off[curr] + next_i;
          const int next_v = adj[e2];

          if (!used[e2]) {
            used[e2] = 1;
            ++used_count;
          }

          prev_v = curr;
          curr = next_v;
          idx_prev = rev_pos[e2];
          ++iterations;
        }

        if (curr == start_u && face.size() >= 3 && iterations < 2 * m_dir) {
          area2 += pts[prev_v].x * pts[start_u].y - pts[start_u].x * pts[prev_v].y;
          const double a = 0.5 * area2;
          // Skip outer face / degenerates.
          if (a > 1e-12) {
            ++dbg_faces_;
            dbg_sum_face_verts_ += static_cast<long long>(face.size());
            triangulate_one_monotone_face(pts, face);
          }
        }
      }
    }

    dbg_used_half_edges_ = used_count;
  }

  void triangulate_one_monotone_face(const std::vector<Point>& pts, const std::vector<int>& face) {
    const int m = static_cast<int>(face.size());
    if (m < 3) return;
    if (m == 3) {
      triangles_.push_back({face[0], face[1], face[2]});
      return;
    }

    // Find top and bottom vertices by (y desc, x asc).
    auto better_top = [&](int a, int b) {
      if (std::abs(pts[a].y - pts[b].y) > detail::kEps) return pts[a].y > pts[b].y;
      return pts[a].x < pts[b].x;
    };
    auto better_bottom = [&](int a, int b) {
      if (std::abs(pts[a].y - pts[b].y) > detail::kEps) return pts[a].y < pts[b].y;
      return pts[a].x < pts[b].x;
    };

    int top = face[0], bottom = face[0];
    for (int v : face) {
      if (better_top(v, top)) top = v;
      if (better_bottom(v, bottom)) bottom = v;
    }

    // Split into two chains from top to bottom along the face cycle.
    const int idx_top = static_cast<int>(std::find(face.begin(), face.end(), top) - face.begin());

    tmp_chain1_.clear();
    tmp_chain2_.clear();
    tmp_chain1_.reserve(m);
    tmp_chain2_.reserve(m);
    // Walk forward (CCW) from top to bottom.
    for (int i = idx_top; ; i = (i + 1) % m) {
      tmp_chain1_.push_back(face[i]);
      if (face[i] == bottom) break;
    }
    // Walk backward (CW) from top to bottom.
    for (int i = idx_top; ; i = (i - 1 + m) % m) {
      tmp_chain2_.push_back(face[i]);
      if (face[i] == bottom) break;
    }

    // Determine which is left vs right by comparing next step from top.
    auto x_next1 = pts[tmp_chain1_.size() > 1 ? tmp_chain1_[1] : bottom].x;
    auto x_next2 = pts[tmp_chain2_.size() > 1 ? tmp_chain2_[1] : bottom].x;
    const bool chain1_is_left = x_next1 < x_next2;

    // Mark chain membership (reuse buffer to avoid per-face allocations).
    for (int v : is_left_touched_) is_left_buf_[v] = 0;
    is_left_touched_.clear();
    const auto& left_chain = chain1_is_left ? tmp_chain1_ : tmp_chain2_;
    const auto& right_chain = chain1_is_left ? tmp_chain2_ : tmp_chain1_;
    for (int v : left_chain) {
      if (!is_left_buf_[v]) {
        is_left_buf_[v] = 1;
        is_left_touched_.push_back(v);
      }
    }

    // Build the sorted-by-y vertex order by merging the two monotone chains
    // (linear time). Both chains include {top, bottom}.
    tmp_sorted_.clear();
    tmp_sorted_.reserve(m);
    tmp_sorted_.push_back(top);
    std::size_t i = 1, j = 1;
    while (i + 1 < left_chain.size() || j + 1 < right_chain.size()) {
      if (i + 1 >= left_chain.size()) {
        tmp_sorted_.push_back(right_chain[j++]);
        continue;
      }
      if (j + 1 >= right_chain.size()) {
        tmp_sorted_.push_back(left_chain[i++]);
        continue;
      }
      if (better_top(left_chain[i], right_chain[j])) {
        tmp_sorted_.push_back(left_chain[i++]);
      } else {
        tmp_sorted_.push_back(right_chain[j++]);
      }
    }
    tmp_sorted_.push_back(bottom);

    // Stack-based monotone triangulation.
    tmp_stack_.clear();
    tmp_stack_.reserve(m);
    tmp_stack_.push_back(tmp_sorted_[0]);
    tmp_stack_.push_back(tmp_sorted_[1]);

    auto same_chain = [&](int a, int b) { return is_left_buf_[a] == is_left_buf_[b]; };

    for (int i = 2; i < m - 1; ++i) {
      const int v = tmp_sorted_[i];
      if (!same_chain(v, tmp_stack_.back())) {
        // Different chains: connect v to all stack vertices.
        while (static_cast<int>(tmp_stack_.size()) > 1) {
          const int u = tmp_stack_.back(); tmp_stack_.pop_back();
          const int w = tmp_stack_.back();
          triangles_.push_back({v, u, w});
        }
        tmp_stack_.pop_back();
        tmp_stack_.push_back(tmp_sorted_[i - 1]);
        tmp_stack_.push_back(v);
      } else {
        // Same chain: pop while diagonal is inside.
        int u = tmp_stack_.back(); tmp_stack_.pop_back();
        while (!tmp_stack_.empty()) {
          const int w = tmp_stack_.back();
          const double c = detail::cross(pts[v], pts[u], pts[w]);
          const bool ok = is_left_buf_[v] ? (c > detail::kEps) : (c < -detail::kEps);
          if (!ok) break;
          triangles_.push_back({v, u, w});
          u = tmp_stack_.back(); tmp_stack_.pop_back();
        }
        tmp_stack_.push_back(u);
        tmp_stack_.push_back(v);
      }
    }

    // Connect bottom to remaining stack.
    const int v = tmp_sorted_[m - 1];  // bottom
    while (static_cast<int>(tmp_stack_.size()) > 1) {
      const int u = tmp_stack_.back(); tmp_stack_.pop_back();
      const int w = tmp_stack_.back();
      triangles_.push_back({v, u, w});
    }
  }
};

}  // namespace reflex_tri


