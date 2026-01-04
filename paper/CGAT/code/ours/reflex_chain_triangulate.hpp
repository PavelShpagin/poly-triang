#pragma once
/**
 * Chain-based polygon triangulation (paper implementation).
 *
 * Target complexity: O(n + r log r), where r is the number of reflex vertices.
 *
 * Notes:
 * - We assume a simple polygon without holes.
 * - We support equal y-coordinates via a deterministic tie-break (y, then x),
 *   but the theoretical paper assumes general position.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <stdexcept>
#include <unordered_set>
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

inline double cross(const Point& a, const Point& b, const Point& c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
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
  std::vector<Triangle> triangulate(std::vector<Point>& pts) {
    triangles_.clear();
    diagonals_.clear();
    diag_set_.clear();
    reflex_count_ = 0;
    dbg_faces_ = 0;
    dbg_sum_face_verts_ = 0;
    dbg_half_edges_ = 0;
    dbg_used_half_edges_ = 0;

    const int n = static_cast<int>(pts.size());
    if (n < 3) return {};
    diag_set_.reserve(static_cast<std::size_t>(2 * n + 8));

    // Ensure indices and CCW orientation.
    for (int i = 0; i < n; ++i) pts[i].index = i;
    if (detail::signed_area(pts) < 0.0) {
      std::reverse(pts.begin(), pts.end());
      for (int i = 0; i < n; ++i) pts[i].index = i;
    }

    // NOTE: The paper assumes general position (distinct y-coordinates).
    // For performance (and to preserve the output-sensitive bound), we do not
    // enforce this via an O(n log n) check here.

    // Count all reflex vertices (paper's r).
    for (int i = 0; i < n; ++i) {
      if (detail::is_reflex_vertex_ccw(pts, i)) ++reflex_count_;
    }

    if (n == 3) {
      triangles_.push_back({0, 1, 2});
      return triangles_;
    }

    if (reflex_count_ == 0) {
      // Convex -> fan triangulation.
      triangles_.reserve(n - 2);
      for (int i = 1; i < n - 1; ++i) triangles_.push_back({0, i, i + 1});
      return triangles_;
    }

    build_vertex_types_and_chains(pts);
    decompose_into_monotone(pts);
    triangulate_monotone_faces(pts);

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

  struct Chain {
    std::vector<int> verts;  // from upper (local max) to lower (local min)
    int curr = 0;            // current edge is verts[curr] -> verts[curr+1]
    int pending = -1;        // pending merge vertex index (or -1)
    int upper = -1;
    int lower = -1;

    void reset_runtime() {
      curr = 0;
      pending = -1;
    }

    void advance(double y, const std::vector<Point>& pts) {
      // Move down the chain until current edge spans y.
      while (curr + 1 < static_cast<int>(verts.size()) &&
             pts[verts[curr + 1]].y > y + detail::kEps) {
        ++curr;
      }
      if (curr + 1 >= static_cast<int>(verts.size())) {
        curr = static_cast<int>(verts.size()) - 2;
      }
      if (curr < 0) curr = 0;
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
  std::vector<std::pair<int, int>> diagonals_;
  // For uniqueness without an O(|D| log |D|) sort pass.
  std::unordered_set<std::uint64_t> diag_set_;
  std::vector<DiagRecord> diag_records_;
  std::vector<PendingRecord> pending_records_;
  std::vector<Triangle> triangles_;
  int reflex_count_ = 0;
  // Debug stats from face extraction (helpful to diagnose degeneracy / bugs).
  int dbg_faces_ = 0;
  long long dbg_sum_face_verts_ = 0;
  std::size_t dbg_half_edges_ = 0;
  std::size_t dbg_used_half_edges_ = 0;

  // Scratch buffers (reused across faces) to avoid per-face allocations in
  // monotone triangulation.
  std::vector<char> is_left_buf_;
  std::vector<int> is_left_touched_;

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
    StatusNode* parent = nullptr;
  };

  std::vector<StatusNode*> status_node_of_chain_;

  static std::uint64_t diag_key(int a, int b) {
    // Undirected key with a<b.
    if (a > b) std::swap(a, b);
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) |
           static_cast<std::uint32_t>(b);
  }

  bool add_diagonal_unique(int u, int v, int event_v, VType event_type, int chosen_chain,
                           int target_v, const char* reason) {
    if (u == v) return false;
    if (u > v) std::swap(u, v);
    const std::uint64_t key = diag_key(u, v);
    if (!diag_set_.insert(key).second) return false;
    diagonals_.push_back({u, v});
    DiagRecord r;
    r.a = u;
    r.b = v;
    r.event_v = event_v;
    r.event_type = static_cast<int>(event_type);
    r.chosen_chain = chosen_chain;
    r.target_v = target_v;
    r.reason = reason;
    diag_records_.push_back(r);
    return true;
  }

  void advance_chain(int cid, const std::vector<Point>& pts) {
    Chain& ch = chains_[cid];
    const int n = static_cast<int>(pts.size());
    // Step the chain pointer down while the next vertex is still above sweep_y_.
    while (ch.curr + 1 < static_cast<int>(ch.verts.size()) &&
           pts[ch.verts[ch.curr + 1]].y > sweep_y_ + detail::kEps) {
      ++ch.curr;
      const int vreg = ch.verts[ch.curr];
      if (types_[vreg] != VType::Regular) continue;
      if (ch.pending == -1) continue;
      // Only the "interior to right" regular case corresponds to helper(e_{i-1})
      // on this chain (PolyPartition's REGULAR/right case):
      // Below(v_i, v_{i-1}).
      const int prev = (vreg - 1 + n) % n;
      if (!detail::below(pts, vreg, prev)) continue;

      add_diagonal_unique(vreg, ch.pending, /*event_v=*/vreg, VType::Regular,
                          /*chosen_chain=*/cid, /*target_v=*/ch.pending,
                          /*reason=*/"regular:pending");
      // After the regular-vertex action, helper(e_i) becomes vreg (not a merge),
      // so we clear the pending merge.
      ch.pending = -1;
    }
    if (ch.curr + 1 >= static_cast<int>(ch.verts.size())) {
      ch.curr = static_cast<int>(ch.verts.size()) - 2;
    }
    if (ch.curr < 0) ch.curr = 0;
  }

  double chain_x_at(int cid, const std::vector<Point>& pts) {
    advance_chain(cid, pts);
    const int a = chains_[cid].verts[chains_[cid].curr];
    const int b = chains_[cid].verts[chains_[cid].curr + 1];
    const auto& p = pts[a];
    const auto& q = pts[b];
    if (std::abs(p.y - q.y) < detail::kEps) return std::min(p.x, q.x);
    const double t = (p.y - sweep_y_) / (p.y - q.y);  // p.y >= sweep_y_ >= q.y (in general position)
    return p.x + t * (q.x - p.x);
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

  void set_parent(StatusNode* child, StatusNode* p) {
    if (child) child->parent = p;
  }

  StatusNode* rotate_merge(StatusNode* a, StatusNode* b) {
    if (!a) return b;
    if (!b) return a;
    if (a->prio < b->prio) {
      a->right = rotate_merge(a->right, b);
      set_parent(a->right, a);
      a->parent = nullptr;
      return a;
    }
    b->left = rotate_merge(a, b->left);
    set_parent(b->left, b);
    b->parent = nullptr;
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
      set_parent(root->right, root);
      a = root;
      a->parent = nullptr;
    } else {
      split_by_chain(root->left, pivot_chain, pts, a, root->left);
      set_parent(root->left, root);
      b = root;
      b->parent = nullptr;
    }
  }

  void treap_insert(StatusNode*& root, StatusNode* node, const std::vector<Point>& pts) {
    if (!root) {
      root = node;
      node->parent = nullptr;
      return;
    }
    if (node->prio < root->prio) {
      split_by_chain(root, node->chain_id, pts, node->left, node->right);
      set_parent(node->left, node);
      set_parent(node->right, node);
      root = node;
      root->parent = nullptr;
      return;
    }
    if (chain_less(node->chain_id, root->chain_id, pts)) {
      treap_insert(root->left, node, pts);
      set_parent(root->left, root);
    } else {
      treap_insert(root->right, node, pts);
      set_parent(root->right, root);
    }
    root->parent = nullptr;
  }

  void treap_erase(StatusNode*& root, int chain_id, const std::vector<Point>& pts) {
    if (!root) return;
    if (root->chain_id == chain_id) {
      StatusNode* old = root;
      root = rotate_merge(root->left, root->right);
      if (root) root->parent = nullptr;
      delete old;
      return;
    }
    if (chain_less(chain_id, root->chain_id, pts)) {
      treap_erase(root->left, chain_id, pts);
      set_parent(root->left, root);
    } else {
      treap_erase(root->right, chain_id, pts);
      set_parent(root->right, root);
    }
    root->parent = nullptr;
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

  void free_status(StatusNode* root) {
    if (!root) return;
    free_status(root->left);
    free_status(root->right);
    delete root;
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
      const double c = detail::cross(pts[p], pts[i], pts[nx]);
      const bool convex = c > detail::kEps;
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
    diag_set_.clear();
    diag_set_.reserve(static_cast<std::size_t>(2 * n + 8));

    // Build extrema event list (size O(r)).
    std::vector<Event> events;
    events.reserve(n);
    for (int i = 0; i < n; ++i) {
      if (types_[i] == VType::Start || types_[i] == VType::Split ||
          types_[i] == VType::End || types_[i] == VType::Merge) {
        events.push_back({i, pts[i].y, pts[i].x, types_[i]});
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

    auto insert_chain = [&](int cid) {
      StatusNode* node = new StatusNode();
      node->chain_id = cid;
      node->prio = next_prio();
      node->left = node->right = node->parent = nullptr;
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
      const double neg_inf = -std::numeric_limits<double>::infinity();
      const double pos_inf = std::numeric_limits<double>::infinity();
      if (ev.type == VType::End || ev.type == VType::Merge) {
        sweep_y_ = std::nextafter(vy, pos_inf);
      } else {
        sweep_y_ = std::nextafter(vy, neg_inf);
      }

      if (ev.type == VType::Start) {
        const int cid = max_to_chain_[v];
        if (cid < 0) throw std::runtime_error("missing left chain for Start");
        insert_chain(cid);
      } else if (ev.type == VType::End) {
        const int rid = min_to_chain_[v];
        if (rid < 0) throw std::runtime_error("missing left chain for End");
        // Ensure all skipped REGULAR vertices on R are processed before we act on pending.
        advance_chain(rid, pts);
        if (chains_[rid].pending != -1) {
          add_diagonal_unique(v, chains_[rid].pending, /*event_v=*/v, ev.type,
                              /*chosen_chain=*/rid, /*target_v=*/chains_[rid].pending,
                              /*reason=*/"end:pending");
          chains_[rid].pending = -1;
        }
        remove_chain(rid);
      } else if (ev.type == VType::Split) {
        const int left_id = predecessor_chain_by_x(status, pts, pts[v].x);
        if (left_id >= 0) {
          // Ensure L is advanced to this sweep height before inspecting pending / slab entry.
          advance_chain(left_id, pts);
          if (chains_[left_id].pending != -1) {
            add_diagonal_unique(v, chains_[left_id].pending, /*event_v=*/v, ev.type,
                                /*chosen_chain=*/left_id, /*target_v=*/chains_[left_id].pending,
                                /*reason=*/"split:pending");
            chains_[left_id].pending = -1;
          } else {
            const int tgt = chain_slab_entry_vertex(left_id, pts);
            add_diagonal_unique(v, tgt, /*event_v=*/v, ev.type,
                                /*chosen_chain=*/left_id, /*target_v=*/tgt,
                                /*reason=*/"split:slab");
          }
        }
        const int cid = max_to_chain_[v];
        if (cid < 0) throw std::runtime_error("missing left chain for Split");
        insert_chain(cid);
      } else if (ev.type == VType::Merge) {
        const int rid = min_to_chain_[v];
        if (rid < 0) throw std::runtime_error("missing left chain for Merge");

        // Ensure all skipped REGULAR vertices on R are processed before we act on pending.
        advance_chain(rid, pts);
        if (chains_[rid].pending != -1) {
          add_diagonal_unique(v, chains_[rid].pending, /*event_v=*/v, ev.type,
                              /*chosen_chain=*/rid, /*target_v=*/chains_[rid].pending,
                              /*reason=*/"merge:term_pending");
          chains_[rid].pending = -1;
        }

        // Remove R from the status structure (edge e_{i-1}).
        remove_chain(rid);

        // After removal, find the chain immediately left of v (edge e_j).
        const int left_id = predecessor_chain_by_x(status, pts, pts[v].x);
        if (left_id >= 0) {
          // Ensure L is advanced to this sweep height before inspecting pending.
          advance_chain(left_id, pts);
          if (chains_[left_id].pending != -1) {
            add_diagonal_unique(v, chains_[left_id].pending, /*event_v=*/v, ev.type,
                                /*chosen_chain=*/left_id, /*target_v=*/chains_[left_id].pending,
                                /*reason=*/"merge:left_pending");
            chains_[left_id].pending = -1;
          }
          chains_[left_id].pending = v;
          pending_records_.push_back({left_id, v, v, "merge:set_pending"});
        }
      } else {
        // Regular vertices are skipped in extrema-only sweep (handled implicitly via advance_chain()).
      }
    }

    free_status(status);
  }

  // --- Phase 3: face extraction + monotone triangulation ---
  void triangulate_monotone_faces(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());

    // Build adjacency (undirected) from polygon boundary + diagonals.
    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; ++i) {
      const int j = (i + 1) % n;
      adj[i].push_back(j);
      adj[j].push_back(i);
    }
    for (const auto& [a, b] : diagonals_) {
      adj[a].push_back(b);
      adj[b].push_back(a);
    }

    // Sort neighbors around each vertex by angle (CCW), without atan2 (faster).
    for (int u = 0; u < n; ++u) {
      auto& nb = adj[u];
      std::sort(nb.begin(), nb.end());
      nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
      if (nb.size() <= 2) continue;  // ordering is irrelevant for degree <= 2
      auto half = [&](int v) -> int {
        const double dx = pts[v].x - pts[u].x;
        const double dy = pts[v].y - pts[u].y;
        if (dy > 0) return 0;
        if (dy < 0) return 1;
        return (dx >= 0) ? 0 : 1;
      };
      std::sort(nb.begin(), nb.end(), [&](int a, int b) {
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
    //
    // Performance note:
    // Avoid unordered_set for half-edge bookkeeping; the total number of directed edges
    // is only O(n + |D|), so a flat visited array indexed by (vertex, neighbor-slot)
    // is significantly faster in practice.
    std::vector<int> off(n + 1, 0);
    for (int u = 0; u < n; ++u) off[u + 1] = off[u] + static_cast<int>(adj[u].size());
    const int m_dir = off[n];
    dbg_half_edges_ = m_dir;

    std::vector<unsigned char> used(m_dir, 0);
    int used_count = 0;

    // Reverse lookup: for directed edge (u -> adj[u][i]), store the index of u in adj[v].
    // This makes face-walk steps O(1) without per-step neighbor scans.
    std::vector<int> rev_pos(m_dir, -1);
    for (int u = 0; u < n; ++u) {
      const auto& nbu = adj[u];
      for (int i = 0; i < static_cast<int>(nbu.size()); ++i) {
        const int v = nbu[i];
        const auto& nbv = adj[v];
        int j = -1;
        for (int t = 0; t < static_cast<int>(nbv.size()); ++t) {
          if (nbv[t] == u) { j = t; break; }
        }
        rev_pos[off[u] + i] = j;
      }
    }

    auto face_area = [&](const std::vector<int>& cyc) -> double {
      double a = 0.0;
      for (int i = 0; i < static_cast<int>(cyc.size()); ++i) {
        const int j = (i + 1) % static_cast<int>(cyc.size());
        a += pts[cyc[i]].x * pts[cyc[j]].y - pts[cyc[j]].x * pts[cyc[i]].y;
      }
      return 0.5 * a;
    };

    std::vector<std::vector<int>> faces_all;
    std::vector<double> face_areas;
    faces_all.reserve(diagonals_.size() + 4);
    face_areas.reserve(diagonals_.size() + 4);

    for (int start_u = 0; start_u < n; ++start_u) {
      const auto& nb0 = adj[start_u];
      for (int start_i = 0; start_i < static_cast<int>(nb0.size()); ++start_i) {
        const int start_e = off[start_u] + start_i;
        if (used[start_e]) continue;

        std::vector<int> face;
        face.reserve(16);
        face.push_back(start_u);

        int prev = start_u;
        int curr = nb0[start_i];

        used[start_e] = 1;
        ++used_count;
        int idx_prev = rev_pos[start_e];

        int iterations = 0;
        while (curr != start_u && iterations < 2 * m_dir) {
          face.push_back(curr);

          const auto& nb = adj[curr];
          if (idx_prev < 0 || idx_prev >= static_cast<int>(nb.size())) break;

          const int next_i = (idx_prev - 1 + static_cast<int>(nb.size())) % static_cast<int>(nb.size());
          const int next_v = nb[next_i];

          const int e2 = off[curr] + next_i;
          if (!used[e2]) {
            used[e2] = 1;
            ++used_count;
          }

          prev = curr;
          curr = next_v;
          idx_prev = rev_pos[e2];
          ++iterations;
        }

        if (curr == start_u && face.size() >= 3 && iterations < 2 * m_dir) {
          const double a = face_area(face);
          faces_all.push_back(std::move(face));
          face_areas.push_back(a);
        }
      }
    }

    dbg_used_half_edges_ = used_count;

    // Identify the outer face robustly as the face with the largest absolute area,
    // then triangulate all remaining (interior) faces after orienting them CCW.
    std::vector<std::vector<int>> faces;
    faces.reserve(faces_all.size());
    if (!faces_all.empty()) {
      std::size_t outer_idx = 0;
      double best_abs = std::abs(face_areas[0]);
      for (std::size_t i = 1; i < faces_all.size(); ++i) {
        const double ab = std::abs(face_areas[i]);
        if (ab > best_abs) {
          best_abs = ab;
          outer_idx = i;
        }
      }
      for (std::size_t i = 0; i < faces_all.size(); ++i) {
        if (i == outer_idx) continue;
        if (std::abs(face_areas[i]) <= 1e-12) continue;
        auto face = std::move(faces_all[i]);
        if (face_areas[i] < 0.0) std::reverse(face.begin(), face.end());
        faces.push_back(std::move(face));
      }
    }

    dbg_faces_ = static_cast<int>(faces.size());
    dbg_sum_face_verts_ = 0;
    for (const auto& f : faces) dbg_sum_face_verts_ += static_cast<long long>(f.size());

    // Triangulate each monotone face.
    if (static_cast<int>(is_left_buf_.size()) != n) {
      is_left_buf_.assign(n, 0);
    } else {
      std::fill(is_left_buf_.begin(), is_left_buf_.end(), 0);
    }
    is_left_touched_.clear();
    for (const auto& f : faces) {
      triangulate_one_monotone_face(pts, f);
    }
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
    const int idx_bot = static_cast<int>(std::find(face.begin(), face.end(), bottom) - face.begin());

    std::vector<int> chain1, chain2;
    // Walk forward (CCW) from top to bottom.
    for (int i = idx_top; ; i = (i + 1) % m) {
      chain1.push_back(face[i]);
      if (face[i] == bottom) break;
    }
    // Walk backward (CW) from top to bottom.
    for (int i = idx_top; ; i = (i - 1 + m) % m) {
      chain2.push_back(face[i]);
      if (face[i] == bottom) break;
    }

    // Determine which is left vs right by comparing next step from top.
    auto x_next1 = pts[chain1.size() > 1 ? chain1[1] : bottom].x;
    auto x_next2 = pts[chain2.size() > 1 ? chain2[1] : bottom].x;
    const bool chain1_is_left = x_next1 < x_next2;

    // Mark chain membership (reuse buffer to avoid per-face allocations).
    for (int v : is_left_touched_) is_left_buf_[v] = 0;
    is_left_touched_.clear();
    const auto& left_chain = chain1_is_left ? chain1 : chain2;
    const auto& right_chain = chain1_is_left ? chain2 : chain1;
    for (int v : left_chain) {
      if (!is_left_buf_[v]) {
        is_left_buf_[v] = 1;
        is_left_touched_.push_back(v);
      }
    }

    // Build the sorted-by-y vertex order by merging the two monotone chains
    // (linear time). Both chains include {top, bottom}.
    std::vector<int> sorted;
    sorted.reserve(m);
    sorted.push_back(top);
    std::size_t i = 1, j = 1;
    while (i + 1 < left_chain.size() || j + 1 < right_chain.size()) {
      if (i + 1 >= left_chain.size()) {
        sorted.push_back(right_chain[j++]);
        continue;
      }
      if (j + 1 >= right_chain.size()) {
        sorted.push_back(left_chain[i++]);
        continue;
      }
      if (better_top(left_chain[i], right_chain[j])) {
        sorted.push_back(left_chain[i++]);
      } else {
        sorted.push_back(right_chain[j++]);
      }
    }
    sorted.push_back(bottom);

    // Stack-based monotone triangulation.
    std::vector<int> st;
    st.push_back(sorted[0]);
    st.push_back(sorted[1]);

    auto same_chain = [&](int a, int b) { return is_left_buf_[a] == is_left_buf_[b]; };

    for (int i = 2; i < m - 1; ++i) {
      const int v = sorted[i];
      if (!same_chain(v, st.back())) {
        // Different chains: connect v to all stack vertices.
        while (static_cast<int>(st.size()) > 1) {
          const int u = st.back(); st.pop_back();
          const int w = st.back();
          triangles_.push_back({v, u, w});
        }
        st.pop_back();
        st.push_back(sorted[i - 1]);
        st.push_back(v);
      } else {
        // Same chain: pop while diagonal is inside.
        int u = st.back(); st.pop_back();
        while (!st.empty()) {
          const int w = st.back();
          const double c = detail::cross(pts[v], pts[u], pts[w]);
          const bool ok = is_left_buf_[v] ? (c > detail::kEps) : (c < -detail::kEps);
          if (!ok) break;
          triangles_.push_back({v, u, w});
          u = st.back(); st.pop_back();
        }
        st.push_back(u);
        st.push_back(v);
      }
    }

    // Connect bottom to remaining stack.
    const int v = sorted[m - 1];  // bottom
    while (static_cast<int>(st.size()) > 1) {
      const int u = st.back(); st.pop_back();
      const int w = st.back();
      triangles_.push_back({v, u, w});
    }
  }
};

}  // namespace reflex_tri


