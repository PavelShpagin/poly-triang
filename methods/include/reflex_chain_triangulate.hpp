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
  // "a below b" w/ deterministic tie-breaker: for equal y, larger x is below.
  if (std::abs(pts[a].y - pts[b].y) > kEps) return pts[a].y < pts[b].y;
  return pts[a].x > pts[b].x;
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
    reflex_count_ = 0;
    dbg_faces_ = 0;
    dbg_sum_face_verts_ = 0;
    dbg_half_edges_ = 0;
    dbg_used_half_edges_ = 0;

    const int n = static_cast<int>(pts.size());
    if (n < 3) return {};

    // Ensure indices and CCW orientation.
    for (int i = 0; i < n; ++i) pts[i].index = i;
    if (detail::signed_area(pts) < 0.0) {
      std::reverse(pts.begin(), pts.end());
      for (int i = 0; i < n; ++i) pts[i].index = i;
    }

    // Paper assumption: general position (no two vertices share the same y).
    // For benchmarks we generate polygons in general position deterministically.
    {
      std::vector<double> ys;
      ys.reserve(n);
      for (const auto& p : pts) ys.push_back(p.y);
      std::sort(ys.begin(), ys.end());
      for (int i = 1; i < n; ++i) {
        if (std::abs(ys[i] - ys[i - 1]) <= detail::kEps) {
          throw std::runtime_error(
              "input violates general position: two vertices have the same y "
              "(paper assumes distinct y-coordinates)");
        }
      }
    }

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

    double x_at(double y, const std::vector<Point>& pts) {
      advance(y, pts);
      const int a = verts[curr];
      const int b = verts[curr + 1];
      const auto& p = pts[a];
      const auto& q = pts[b];
      if (std::abs(p.y - q.y) < detail::kEps) return std::min(p.x, q.x);
      const double t = (p.y - y) / (p.y - q.y);  // p.y >= y >= q.y
      return p.x + t * (q.x - p.x);
    }
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
  std::vector<DiagRecord> diag_records_;
  std::vector<PendingRecord> pending_records_;
  std::vector<Triangle> triangles_;
  int reflex_count_ = 0;
  // Debug stats from face extraction (helpful to diagnose degeneracy / bugs).
  int dbg_faces_ = 0;
  long long dbg_sum_face_verts_ = 0;
  std::size_t dbg_half_edges_ = 0;
  std::size_t dbg_used_half_edges_ = 0;

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

  bool chain_less(int a, int b, const std::vector<Point>& pts) {
    if (a == b) return false;
    const double xa = chains_[a].x_at(sweep_y_, pts);
    const double xb = chains_[b].x_at(sweep_y_, pts);
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
      const double xcur = chains_[root->chain_id].x_at(sweep_y_, pts);
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
    const double x = chains_[cid].x_at(sweep_y_, pts);
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

    // --- Regular-vertex handling (implicit) ---
    //
    // The chain/event sweep skips explicit REGULAR events, but we still must
    // replicate the classical behavior that resolves a pending MERGE helper on
    // an edge as soon as we pass the next REGULAR vertex on the same chain
    // whose interior lies to the right (PolyPartition's REGULAR/right case).
    //
    // Without this, a MERGE can remain "pending" too long and a later extremum
    // may connect to it across boundary edges, producing invalid diagonals.
    auto advance_chain_with_regular_resolution = [&](int cid, double y) {
      Chain& ch = chains_[cid];
      // Step the chain pointer down while the next vertex is still above y.
      while (ch.curr + 1 < static_cast<int>(ch.verts.size()) &&
             pts[ch.verts[ch.curr + 1]].y > y + detail::kEps) {
        ++ch.curr;
        const int vreg = ch.verts[ch.curr];
        if (types_[vreg] != VType::Regular) continue;
        if (ch.pending == -1) continue;
        // Only the "interior to right" regular case corresponds to helper(e_{i-1})
        // on this chain. This is exactly the condition used in PolyPartition:
        // Below(v_i, v_{i-1}).
        const int prev = (vreg - 1 + n) % n;
        if (!detail::below(pts, vreg, prev)) continue;

        diagonals_.push_back({vreg, ch.pending});
        DiagRecord r;
        r.a = std::min(vreg, ch.pending);
        r.b = std::max(vreg, ch.pending);
        r.event_v = vreg;
        r.event_type = static_cast<int>(VType::Regular);
        r.chosen_chain = cid;
        r.target_v = ch.pending;
        r.reason = "regular:pending";
        diag_records_.push_back(r);

        // After the regular-vertex action, helper(e_i) becomes vreg (not a merge),
        // so we clear the pending merge.
        ch.pending = -1;
      }
      if (ch.curr + 1 >= static_cast<int>(ch.verts.size())) {
        ch.curr = static_cast<int>(ch.verts.size()) - 2;
      }
      if (ch.curr < 0) ch.curr = 0;
    };

    // Build extrema event list.
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

    // Status structure:
    // For correctness, we follow the reference implementation logic:
    // maintain a simple list of active left-chains and find the chain left of a vertex
    // by scanning for the maximal x_at(sweep_y_) strictly less than v.x.
    //
    // This avoids relying on treap ordering when chain pointers advance lazily.
    std::vector<int> active;
    active.reserve(chains_.size());

    for (const auto& ev : events) {
      const int v = ev.v;
      // Sweep convention (paper): use y-eps for maxima, y+eps for minima.
      if (ev.type == VType::Start || ev.type == VType::Split) sweep_y_ = pts[v].y - 1e-9;
      else sweep_y_ = pts[v].y + 1e-9;

      // IMPORTANT: advance all active chains to this sweep height once. This is
      // where we implicitly process any skipped REGULAR vertices and may
      // generate diagonals that resolve pending merges early.
      for (int cid : active) {
        advance_chain_with_regular_resolution(cid, sweep_y_);
      }

      // Build sorted status order by x at sweep height.
      auto sorted_active = [&]() {
        std::vector<std::pair<double, int>> xs;
        xs.reserve(active.size());
        for (int cid : active) xs.push_back({chains_[cid].x_at(sweep_y_, pts), cid});
        std::sort(xs.begin(), xs.end(), [](const auto& a, const auto& b) {
          if (a.first < b.first) return true;
          if (b.first < a.first) return false;
          return a.second < b.second;
        });
        return xs;
      };

      auto predecessor_of_x = [&](double x_query) -> int {
        const auto xs = sorted_active();
        int best = -1;
        for (const auto& [x, cid] : xs) {
          if (x < x_query) best = cid;
          else break;
        }
        return best;
      };

      auto predecessor_of_chain = [&](int chain_id) -> int {
        const auto xs = sorted_active();
        for (size_t i = 0; i < xs.size(); ++i) {
          if (xs[i].second == chain_id) {
            return (i == 0) ? -1 : xs[i - 1].second;
          }
        }
        return -1;
      };

      auto insert_chain = [&](int cid) {
        chains_[cid].reset_runtime();
        active.push_back(cid);
      };

      auto remove_chain = [&](int cid) {
        for (size_t i = 0; i < active.size(); ++i) {
          if (active[i] == cid) {
            active[i] = active.back();
            active.pop_back();
            return;
          }
        }
      };

      if (ev.type == VType::Start) {
        const int cid = max_to_chain_[v];
        if (cid < 0) throw std::runtime_error("missing left chain for Start");
        insert_chain(cid);
      } else if (ev.type == VType::End) {
        const int rid = min_to_chain_[v];
        if (rid < 0) throw std::runtime_error("missing left chain for End");
        if (chains_[rid].pending != -1) {
          diagonals_.push_back({v, chains_[rid].pending});
          DiagRecord r;
          r.a = std::min(v, chains_[rid].pending);
          r.b = std::max(v, chains_[rid].pending);
          r.event_v = v;
          r.event_type = static_cast<int>(ev.type);
          r.chosen_chain = rid;
          r.target_v = chains_[rid].pending;
          r.reason = "end:pending";
          diag_records_.push_back(r);
          chains_[rid].pending = -1;
        }
        remove_chain(rid);
      } else if (ev.type == VType::Split) {
        const int left_id = predecessor_of_x(pts[v].x);
        if (left_id >= 0) {
          if (chains_[left_id].pending != -1) {
            diagonals_.push_back({v, chains_[left_id].pending});
            DiagRecord r;
            r.a = std::min(v, chains_[left_id].pending);
            r.b = std::max(v, chains_[left_id].pending);
            r.event_v = v;
            r.event_type = static_cast<int>(ev.type);
            r.chosen_chain = left_id;
            r.target_v = chains_[left_id].pending;
            r.reason = "split:pending";
            diag_records_.push_back(r);
            chains_[left_id].pending = -1;
          } else {
            diagonals_.push_back({v, chains_[left_id].slab_entry_vertex()});
            DiagRecord r;
            const int tgt = chains_[left_id].slab_entry_vertex();
            r.a = std::min(v, tgt);
            r.b = std::max(v, tgt);
            r.event_v = v;
            r.event_type = static_cast<int>(ev.type);
            r.chosen_chain = left_id;
            r.target_v = tgt;
            r.reason = "split:slab";
            diag_records_.push_back(r);
          }
        }

        const int cid = max_to_chain_[v];
        if (cid < 0) throw std::runtime_error("missing left chain for Split");
        insert_chain(cid);
      } else if (ev.type == VType::Merge) {
        const int rid = min_to_chain_[v];
        if (rid < 0) throw std::runtime_error("missing left chain for Merge");

        if (chains_[rid].pending != -1) {
          diagonals_.push_back({v, chains_[rid].pending});
          DiagRecord r;
          r.a = std::min(v, chains_[rid].pending);
          r.b = std::max(v, chains_[rid].pending);
          r.event_v = v;
          r.event_type = static_cast<int>(ev.type);
          r.chosen_chain = rid;
          r.target_v = chains_[rid].pending;
          r.reason = "merge:term_pending";
          diag_records_.push_back(r);
          chains_[rid].pending = -1;
        }
        // Delete R from status before searching the chain directly left of v.
        remove_chain(rid);

        // Match PolyPartition: search by x(v) (with R removed) to find the chain
        // immediately left of v at this sweep height.
        const int left_id = predecessor_of_x(pts[v].x);
        if (left_id >= 0) {
          if (chains_[left_id].pending != -1) {
            diagonals_.push_back({v, chains_[left_id].pending});
            DiagRecord r;
            r.a = std::min(v, chains_[left_id].pending);
            r.b = std::max(v, chains_[left_id].pending);
            r.event_v = v;
            r.event_type = static_cast<int>(ev.type);
            r.chosen_chain = left_id;
            r.target_v = chains_[left_id].pending;
            r.reason = "merge:left_pending";
            diag_records_.push_back(r);
            chains_[left_id].pending = -1;
          }
          chains_[left_id].pending = v;
          pending_records_.push_back({left_id, v, v, "merge:set_pending"});
        }
      } else {
        // Regular vertices are skipped in extrema-only sweep (handled implicitly by chain pointers).
      }
    }

    // Deduplicate diagonals (can occur if events coincide in degenerate inputs).
    for (auto& d : diagonals_) {
      if (d.first > d.second) std::swap(d.first, d.second);
    }
    std::sort(diagonals_.begin(), diagonals_.end());
    diagonals_.erase(std::unique(diagonals_.begin(), diagonals_.end()), diagonals_.end());
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

    // Sort neighbors around each vertex by angle (CCW).
    for (int u = 0; u < n; ++u) {
      auto& nb = adj[u];
      std::sort(nb.begin(), nb.end());
      nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
      std::sort(nb.begin(), nb.end(), [&](int a, int b) {
        const double aa = std::atan2(pts[a].y - pts[u].y, pts[a].x - pts[u].x);
        const double ab = std::atan2(pts[b].y - pts[u].y, pts[b].x - pts[u].x);
        return aa < ab;
      });
    }

    // Face extraction: traverse directed half-edges in the planar embedding.
    // For directed edge (prev -> curr), the next vertex is the clockwise predecessor
    // of `prev` in `curr`'s CCW-sorted neighbor list. This follows the face on the
    // left of the directed edge and yields CCW-oriented interior faces.
    struct PairHash {
      std::size_t operator()(const std::pair<int, int>& p) const noexcept {
        return (static_cast<std::size_t>(p.first) << 32) ^ static_cast<std::size_t>(p.second);
      }
    };

    std::unordered_set<std::pair<int, int>, PairHash> half_edges;
    half_edges.reserve(static_cast<std::size_t>(2 * (n + diagonals_.size()) + 8));
    for (int i = 0; i < n; ++i) {
      const int j = (i + 1) % n;
      half_edges.insert({i, j});
      half_edges.insert({j, i});
    }
    for (const auto& [a, b] : diagonals_) {
      half_edges.insert({a, b});
      half_edges.insert({b, a});
    }
    dbg_half_edges_ = half_edges.size();

    std::unordered_set<std::pair<int, int>, PairHash> used;
    used.reserve(half_edges.size() * 2);

    auto face_area = [&](const std::vector<int>& cyc) -> double {
      double a = 0.0;
      for (int i = 0; i < static_cast<int>(cyc.size()); ++i) {
        const int j = (i + 1) % static_cast<int>(cyc.size());
        a += pts[cyc[i]].x * pts[cyc[j]].y - pts[cyc[j]].x * pts[cyc[i]].y;
      }
      return 0.5 * a;
    };

    std::vector<std::vector<int>> faces;
    faces.reserve(diagonals_.size() + 4);

    for (const auto& start_edge : half_edges) {
      if (used.find(start_edge) != used.end()) continue;
      const int start_u = start_edge.first;
      const int start_v = start_edge.second;

      std::vector<int> face;
      face.push_back(start_u);

      int prev = start_u;
      int curr = start_v;
      used.insert({prev, curr});

      int iterations = 0;
      while (curr != start_u && iterations < 2 * n) {
        face.push_back(curr);
        const auto& neighbors = adj[curr];
        int idx = -1;
        for (int i = 0; i < static_cast<int>(neighbors.size()); ++i) {
          if (neighbors[i] == prev) { idx = i; break; }
        }
        if (idx < 0) break;
        const int next_v =
            neighbors[(idx - 1 + static_cast<int>(neighbors.size())) % static_cast<int>(neighbors.size())];
        used.insert({curr, next_v});
        prev = curr;
        curr = next_v;
        ++iterations;
      }

      if (face.size() >= 3 && iterations < 2 * n) {
        if (face_area(face) > 1e-9) {
          faces.push_back(std::move(face));
        }
      }
    }
    dbg_used_half_edges_ = used.size();
    dbg_faces_ = static_cast<int>(faces.size());
    dbg_sum_face_verts_ = 0;
    for (const auto& f : faces) dbg_sum_face_verts_ += static_cast<long long>(f.size());

    // Triangulate each monotone face.
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

    std::vector<char> is_left(pts.size(), 0);
    const auto& left_chain = chain1_is_left ? chain1 : chain2;
    const auto& right_chain = chain1_is_left ? chain2 : chain1;
    for (int v : left_chain) is_left[v] = 1;
    for (int v : right_chain) is_left[v] = 0;

    // Sort face vertices by decreasing y (then x asc).
    std::vector<int> sorted = face;
    std::sort(sorted.begin(), sorted.end(), [&](int a, int b) { return better_top(a, b); });

    // Stack-based monotone triangulation.
    std::vector<int> st;
    st.push_back(sorted[0]);
    st.push_back(sorted[1]);

    auto same_chain = [&](int a, int b) { return is_left[a] == is_left[b]; };

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
          const bool ok = is_left[v] ? (c > detail::kEps) : (c < -detail::kEps);
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


