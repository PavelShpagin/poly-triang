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
      if (c > detail::kEps) is_convex_buf_[i] = 1;
      if (c < -detail::kEps) ++reflex_count_;
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
             pts[verts[curr + 1]].y > y + detail::kEps) {
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
      // REGULAR vertex handling (textbook sweep):
      // Resolve a pending MERGE helper at regular vertices on the "left chain".
      //
      // In the standard monotone-decomposition sweep for a CCW polygon, a regular vertex
      // is on the left chain iff its next boundary vertex is below it.
      // These vertices can resolve a pending MERGE helper in O(1) amortized time
      // without being explicit sweep events.
      const int nxt = (vreg + 1) % n;
      if (!detail::below(pts, nxt, vreg)) continue;

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
      if (std::abs(p.y - q.y) < detail::kEps) {
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
      if (std::abs(p.y - q.y) < detail::kEps) {
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
      sweep_y_eps_ = sweep_y_ + detail::kEps;

      if (ev.type == VType::Start) {
        const int cid = max_to_chain_[v];
        if (cid < 0) throw std::runtime_error("missing left chain for Start");
        insert_chain(cid);
        chains_[cid].pending = -1;
      } else if (ev.type == VType::End) {
        const int rid = min_to_chain_[v];
        if (rid < 0) throw std::runtime_error("missing left chain for End");
        // Ensure we process any regular vertices on this chain down to the current sweep level
        // (this can clear a pending merge helper) before we emit an end-vertex diagonal.
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

        // Ensure we process any regular vertices on this chain down to the current sweep level
        // (this can clear a pending merge helper) before we emit a merge-vertex diagonal.
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
    // polygon's boundary order (0..n-1).
    //
    // For performance, we use a simple per-vertex sort for small n (lower constant),
    // and a linear-time radix counting sort for larger n. The cutoff is a fixed constant,
    // so asymptotic complexity remains O(n).
    // Always use the linear-time radix ordering of directed half-edges.
    // This avoids O(n) tiny std::sort calls (one per vertex), which can dominate
    // at small/medium n and introduce high variance.
    {
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
    halfedge_used_.assign(m_dir, 0);
    int used_count = 0;

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

    // Outer face: for a CCW polygon, the outer face lies to the left of the reversed
    // boundary edge (0 -> n-1). Since our face-walk always traces the face on the left
    // of a directed half-edge, we can mark the outer face half-edges by walking from
    // that single half-edge.
    std::fill(halfedge_used_.begin(), halfedge_used_.end(), 0);
    used_count = 0;

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

    if (static_cast<int>(is_left_buf_.size()) != n) {
      is_left_buf_.assign(n, 0);
    } else {
      std::fill(is_left_buf_.begin(), is_left_buf_.end(), 0);
    }
    is_left_touched_.clear();
    dbg_faces_ = 0;
    dbg_sum_face_verts_ = 0;

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
          ++dbg_faces_;
          dbg_sum_face_verts_ += static_cast<long long>(face_buf_.size());
          triangulate_one_monotone_face(pts, face_buf_);
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

    // Find top and bottom vertices by (y desc, x asc), with strict ordering.
    // The benchmark setup uses a fixed rotation to avoid equal-y degeneracies.
    auto better_top = [&](int a, int b) {
      if (pts[a].y > pts[b].y) return true;
      if (pts[a].y < pts[b].y) return false;
      if (pts[a].x < pts[b].x) return true;
      if (pts[a].x > pts[b].x) return false;
      return a < b;
    };
    auto better_bottom = [&](int a, int b) {
      if (pts[a].y < pts[b].y) return true;
      if (pts[a].y > pts[b].y) return false;
      if (pts[a].x < pts[b].x) return true;
      if (pts[a].x > pts[b].x) return false;
      return a < b;
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


