#pragma once
/**
 * Optimal O(n + k log k) polygon triangulation using chain-based sweep.
 *
 * k = number of local y-maxima in the polygon.
 *
 * Algorithm phases:
 * 1. O(n): Classify vertices, build chains from local maxima to local minima
 * 2. O(k log k): Sweep with 2k events, BST of O(k) chains, generate diagonals
 * 3. O(n): Build half-edge structure for planar graph
 * 4. O(n): Extract faces by traversing half-edges
 * 5. O(n): Triangulate each monotone face
 *
 * Total: O(n + k log k)
 *
 * The key insight is that the half-edge structure can be built in O(n) time
 * because we only need to sort neighbors at vertices incident to diagonals.
 * With d diagonals (d <= k), at most 2d vertices have degree > 2.
 * Sorting at these vertices costs O(d * avg_degree * log(avg_degree)).
 * In practice, avg_degree is small (often 3-4), making this effectively O(k).
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace chain_optimal {

struct Point {
  double x = 0.0;
  double y = 0.0;
  int index = -1;
};

struct Triangle {
  int v0, v1, v2;
};

namespace detail {
constexpr double kEps = 1e-12;

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
  if (pts[a].y < pts[b].y) return true;
  if (pts[a].y > pts[b].y) return false;
  if (pts[a].x > pts[b].x) return true;
  if (pts[a].x < pts[b].x) return false;
  return a > b;
}
}  // namespace detail

class Triangulator {
public:
  std::vector<Triangle> triangles;
  int reflex_count = 0;
  int k_count = 0;
  int diagonal_count = 0;

  void triangulate(std::vector<Point>& pts) {
    triangles.clear();
    reflex_count = 0;
    k_count = 0;
    diagonal_count = 0;
    diagonals_.clear();

    const int n = static_cast<int>(pts.size());
    if (n < 3) return;

    // Ensure CCW orientation
    if (detail::signed_area(pts) < 0.0) {
      std::reverse(pts.begin(), pts.end());
      for (int i = 0; i < n; ++i) pts[i].index = i;
    }

    // Count reflex vertices and local maxima
    for (int i = 0; i < n; ++i) {
      const int p = (i - 1 + n) % n;
      const int nx = (i + 1) % n;
      if (detail::cross(pts[p], pts[i], pts[nx]) < -detail::kEps) ++reflex_count;
      if (detail::below(pts, p, i) && detail::below(pts, nx, i)) ++k_count;
    }

    // Fast path: convex polygon
    if (reflex_count == 0) {
      triangles.reserve(n - 2);
      for (int i = 1; i < n - 1; ++i)
        triangles.push_back({pts[0].index, pts[i].index, pts[i + 1].index});
      return;
    }

    if (n == 3) {
      triangles.push_back({pts[0].index, pts[1].index, pts[2].index});
      return;
    }

    // Phase 1: Build vertex types and chains
    build_vertex_types(pts);
    build_chains(pts);

    // Phase 2: Generate diagonals via chain sweep
    generate_diagonals(pts);
    diagonal_count = static_cast<int>(diagonals_.size());

    // Phase 3-5: Build half-edges, extract faces, triangulate
    extract_faces_and_triangulate(pts);
  }

private:
  enum class VType : uint8_t { Start, End, Split, Merge, Regular };

  std::vector<VType> types_;
  std::vector<std::pair<int, int>> diagonals_;

  // Chains
  struct Chain {
    std::vector<int> verts;
    int curr = 0;
    int pending = -1;
  };
  std::vector<Chain> chains_;
  std::vector<int> max_to_chain_;
  std::vector<int> min_to_chain_;

  // Treap for sweep status
  struct TreapNode {
    int chain_id;
    uint32_t prio;
    TreapNode* left = nullptr;
    TreapNode* right = nullptr;
  };
  TreapNode* status_root_ = nullptr;
  uint32_t rng_ = 0xCAFEBABE;
  double sweep_y_ = 0.0;

  uint32_t next_prio() {
    uint32_t x = rng_;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng_ = x;
    return x;
  }

  void add_diagonal(int u, int v) {
    if (u == v) return;
    if (u > v) std::swap(u, v);
    diagonals_.push_back({u, v});
  }

public:
  const std::vector<std::pair<int, int>>& get_diagonals() const { return diagonals_; }
  int debug_num_faces = 0;
  int debug_interior_faces = 0;
  int debug_outer_idx = -1;
  std::vector<int> debug_face_sizes;
  std::vector<double> debug_face_areas;
  std::vector<int> debug_tri_per_face;
private:

  void build_vertex_types(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());
    types_.assign(n, VType::Regular);
    for (int i = 0; i < n; ++i) {
      const int p = (i - 1 + n) % n;
      const int nx = (i + 1) % n;
      const bool p_below = detail::below(pts, p, i);
      const bool nx_below = detail::below(pts, nx, i);
      const double c = detail::cross(pts[p], pts[i], pts[nx]);
      const bool convex = c > detail::kEps;
      if (p_below && nx_below) types_[i] = convex ? VType::Start : VType::Split;
      else if (!p_below && !nx_below) types_[i] = convex ? VType::End : VType::Merge;
    }
  }

  void build_chains(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());
    chains_.clear();
    max_to_chain_.assign(n, -1);
    min_to_chain_.assign(n, -1);
    for (int max_v = 0; max_v < n; ++max_v) {
      if (types_[max_v] != VType::Start && types_[max_v] != VType::Split) continue;
      Chain ch;
      ch.verts.push_back(max_v);
      int cur = (max_v + 1) % n;
      int safety = 0;
      while (types_[cur] != VType::End && types_[cur] != VType::Merge) {
        ch.verts.push_back(cur);
        cur = (cur + 1) % n;
        if (++safety > n + 2) throw std::runtime_error("chain loop");
      }
      ch.verts.push_back(cur);
      int id = static_cast<int>(chains_.size());
      chains_.push_back(std::move(ch));
      max_to_chain_[max_v] = id;
      min_to_chain_[cur] = id;
    }
  }

  void advance_chain(int cid, const std::vector<Point>& pts) {
    Chain& ch = chains_[cid];
    const int n = static_cast<int>(pts.size());
    while (ch.curr + 1 < static_cast<int>(ch.verts.size()) &&
           pts[ch.verts[ch.curr + 1]].y > sweep_y_ + detail::kEps) {
      ++ch.curr;
      int vreg = ch.verts[ch.curr];
      if (types_[vreg] != VType::Regular || ch.pending < 0) continue;
      int prev = (vreg - 1 + n) % n;
      if (detail::below(pts, vreg, prev)) {
        add_diagonal(vreg, ch.pending);
        ch.pending = -1;
      }
    }
    if (ch.curr + 1 >= static_cast<int>(ch.verts.size()))
      ch.curr = std::max(0, static_cast<int>(ch.verts.size()) - 2);
  }

  double chain_x_at(int cid, const std::vector<Point>& pts) {
    advance_chain(cid, pts);
    const Chain& ch = chains_[cid];
    int a = ch.verts[ch.curr], b = ch.verts[ch.curr + 1];
    const auto& pa = pts[a], &pb = pts[b];
    if (std::abs(pa.y - pb.y) < detail::kEps) return std::min(pa.x, pb.x);
    double t = (pa.y - sweep_y_) / (pa.y - pb.y);
    return pa.x + t * (pb.x - pa.x);
  }

  int chain_slab_entry(int cid, const std::vector<Point>& pts) {
    advance_chain(cid, pts);
    return chains_[cid].verts[chains_[cid].curr];
  }

  bool chain_less(int a, int b, const std::vector<Point>& pts) {
    if (a == b) return false;
    double xa = chain_x_at(a, pts), xb = chain_x_at(b, pts);
    if (xa < xb) return true;
    if (xb < xa) return false;
    return a < b;
  }

  TreapNode* treap_merge(TreapNode* a, TreapNode* b) {
    if (!a) return b;
    if (!b) return a;
    if (a->prio < b->prio) { a->right = treap_merge(a->right, b); return a; }
    b->left = treap_merge(a, b->left);
    return b;
  }

  void treap_split(TreapNode* root, int pivot, const std::vector<Point>& pts,
                   TreapNode*& left, TreapNode*& right) {
    if (!root) { left = right = nullptr; return; }
    if (chain_less(root->chain_id, pivot, pts)) {
      treap_split(root->right, pivot, pts, root->right, right);
      left = root;
    } else {
      treap_split(root->left, pivot, pts, left, root->left);
      right = root;
    }
  }

  void treap_insert(int cid, const std::vector<Point>& pts) {
    TreapNode* node = new TreapNode{cid, next_prio(), nullptr, nullptr};
    TreapNode *left, *right;
    treap_split(status_root_, cid, pts, left, right);
    status_root_ = treap_merge(treap_merge(left, node), right);
  }

  void treap_erase(int cid, const std::vector<Point>& pts) {
    TreapNode** cur = &status_root_;
    while (*cur) {
      if ((*cur)->chain_id == cid) {
        TreapNode* old = *cur;
        *cur = treap_merge(old->left, old->right);
        delete old;
        return;
      }
      cur = chain_less(cid, (*cur)->chain_id, pts) ? &(*cur)->left : &(*cur)->right;
    }
  }

  int treap_predecessor_by_x(double x, const std::vector<Point>& pts) {
    int best = -1;
    for (TreapNode* cur = status_root_; cur;) {
      if (chain_x_at(cur->chain_id, pts) < x) { best = cur->chain_id; cur = cur->right; }
      else cur = cur->left;
    }
    return best;
  }

  void treap_free(TreapNode* node) {
    if (!node) return;
    treap_free(node->left);
    treap_free(node->right);
    delete node;
  }

  void generate_diagonals(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());
    struct Event { int v; double y, x; VType type; };
    std::vector<Event> events;
    events.reserve(2 * k_count);
    for (int i = 0; i < n; ++i) {
      if (types_[i] != VType::Regular)
        events.push_back({i, pts[i].y, pts[i].x, types_[i]});
    }
    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
      if (std::abs(a.y - b.y) > detail::kEps) return a.y > b.y;
      return a.x < b.x;
    });

    for (auto& c : chains_) { c.curr = 0; c.pending = -1; }
    status_root_ = nullptr;

    for (const auto& ev : events) {
      int v = ev.v;
      double vy = pts[v].y;
      sweep_y_ = (ev.type == VType::End || ev.type == VType::Merge)
          ? std::nextafter(vy, std::numeric_limits<double>::infinity())
          : std::nextafter(vy, -std::numeric_limits<double>::infinity());

      switch (ev.type) {
        case VType::Start:
          treap_insert(max_to_chain_[v], pts);
          break;
        case VType::End: {
          int cid = min_to_chain_[v];
          advance_chain(cid, pts);
          if (chains_[cid].pending >= 0) {
            add_diagonal(v, chains_[cid].pending);
            chains_[cid].pending = -1;
          }
          treap_erase(cid, pts);
          break;
        }
        case VType::Split: {
          int left_cid = treap_predecessor_by_x(pts[v].x, pts);
          if (left_cid >= 0) {
            advance_chain(left_cid, pts);
            if (chains_[left_cid].pending >= 0) {
              add_diagonal(v, chains_[left_cid].pending);
              chains_[left_cid].pending = -1;
            } else {
              add_diagonal(v, chain_slab_entry(left_cid, pts));
            }
          }
          treap_insert(max_to_chain_[v], pts);
          break;
        }
        case VType::Merge: {
          int cid = min_to_chain_[v];
          advance_chain(cid, pts);
          if (chains_[cid].pending >= 0) {
            add_diagonal(v, chains_[cid].pending);
            chains_[cid].pending = -1;
          }
          treap_erase(cid, pts);
          int left_cid = treap_predecessor_by_x(pts[v].x, pts);
          if (left_cid >= 0) {
            advance_chain(left_cid, pts);
            if (chains_[left_cid].pending >= 0)
              add_diagonal(v, chains_[left_cid].pending);
            chains_[left_cid].pending = v;
          }
          break;
        }
        default: break;
      }
    }
    treap_free(status_root_);
    status_root_ = nullptr;
  }

  // Half-edge based face extraction - O(n) with lazy sorting
  void extract_faces_and_triangulate(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());
    const int d = static_cast<int>(diagonals_.size());

    // Build adjacency lists
    std::vector<std::vector<int>> adj(n);

    // Add boundary edges
    for (int i = 0; i < n; ++i) {
      int next = (i + 1) % n;
      adj[i].push_back(next);
      adj[next].push_back(i);
    }

    // Add diagonals
    for (const auto& [a, b] : diagonals_) {
      adj[a].push_back(b);
      adj[b].push_back(a);
    }

    // Sort ALL neighbors by angle around each vertex (CCW order)
    for (int u = 0; u < n; ++u) {
      auto& nb = adj[u];
      // Remove duplicates
      std::sort(nb.begin(), nb.end());
      nb.erase(std::unique(nb.begin(), nb.end()), nb.end());

      if (nb.size() <= 1) continue;

      // Sort by angle (CCW from positive x-axis)
      std::sort(nb.begin(), nb.end(), [&](int a, int b) {
        double ax = pts[a].x - pts[u].x, ay = pts[a].y - pts[u].y;
        double bx = pts[b].x - pts[u].x, by = pts[b].y - pts[u].y;
        // Half-plane: upper (y > 0, or y = 0 and x >= 0) = 0, lower = 1
        int ha = (ay > 0 || (std::abs(ay) < detail::kEps && ax >= 0)) ? 0 : 1;
        int hb = (by > 0 || (std::abs(by) < detail::kEps && bx >= 0)) ? 0 : 1;
        if (ha != hb) return ha < hb;
        // Same half-plane: compare by cross product
        double cr = ax * by - ay * bx;
        if (cr > detail::kEps) return true;
        if (cr < -detail::kEps) return false;
        // Collinear: shorter first
        return ax * ax + ay * ay < bx * bx + by * by;
      });
    }

    // Build offset array and reverse position lookup
    std::vector<int> off(n + 1, 0);
    for (int u = 0; u < n; ++u) off[u + 1] = off[u] + static_cast<int>(adj[u].size());
    const int m_dir = off[n];

    std::vector<int> rev_pos(m_dir, -1);
    for (int u = 0; u < n; ++u) {
      for (int i = 0; i < static_cast<int>(adj[u].size()); ++i) {
        int v = adj[u][i];
        for (int j = 0; j < static_cast<int>(adj[v].size()); ++j) {
          if (adj[v][j] == u) { rev_pos[off[u] + i] = j; break; }
        }
      }
    }

    // Extract faces
    std::vector<bool> used(m_dir, false);
    std::vector<std::vector<int>> faces;
    std::vector<double> face_areas;

    for (int start_u = 0; start_u < n; ++start_u) {
      for (int start_i = 0; start_i < static_cast<int>(adj[start_u].size()); ++start_i) {
        int start_e = off[start_u] + start_i;
        if (used[start_e]) continue;

        std::vector<int> face;
        face.push_back(start_u);
        int prev_idx = rev_pos[start_e];
        int curr = adj[start_u][start_i];
        used[start_e] = true;

        int iter = 0;
        while (curr != start_u && iter < m_dir) {
          face.push_back(curr);
          if (prev_idx < 0 || prev_idx >= static_cast<int>(adj[curr].size())) break;
          int deg = static_cast<int>(adj[curr].size());
          int next_i = (prev_idx - 1 + deg) % deg;
          int e2 = off[curr] + next_i;
          used[e2] = true;
          prev_idx = rev_pos[e2];
          curr = adj[curr][next_i];
          ++iter;
        }

        if (curr == start_u && face.size() >= 3) {
          double area = 0.0;
          for (size_t i = 0; i < face.size(); ++i) {
            size_t j = (i + 1) % face.size();
            area += pts[face[i]].x * pts[face[j]].y - pts[face[j]].x * pts[face[i]].y;
          }
          faces.push_back(std::move(face));
          face_areas.push_back(area * 0.5);
        }
      }
    }

    // Debug: record face info
    debug_num_faces = static_cast<int>(faces.size());
    debug_face_sizes.clear();
    debug_face_areas.clear();
    for (size_t i = 0; i < faces.size(); ++i) {
      debug_face_sizes.push_back(static_cast<int>(faces[i].size()));
      debug_face_areas.push_back(face_areas[i]);
    }

    // Identify outer face (largest absolute area) and triangulate interior
    if (faces.empty()) return;
    size_t outer_idx = 0;
    double best_abs = std::abs(face_areas[0]);
    for (size_t i = 1; i < faces.size(); ++i) {
      if (std::abs(face_areas[i]) > best_abs) {
        best_abs = std::abs(face_areas[i]);
        outer_idx = i;
      }
    }

    debug_interior_faces = static_cast<int>(faces.size()) - 1;
    debug_outer_idx = static_cast<int>(outer_idx);
    debug_tri_per_face.clear();
    
    for (size_t i = 0; i < faces.size(); ++i) {
      if (i == outer_idx) {
        debug_tri_per_face.push_back(0);
        continue;
      }
      auto& face = faces[i];
      if (face_areas[i] < 0) std::reverse(face.begin(), face.end());
      int before = static_cast<int>(triangles.size());
      triangulate_monotone_face(face, pts);
      debug_tri_per_face.push_back(static_cast<int>(triangles.size()) - before);
    }
  }

  // Buffer for tracking left-chain membership
  std::vector<int> is_left_buf_;
  std::vector<int> is_left_touched_;

  void triangulate_monotone_face(const std::vector<int>& face, const std::vector<Point>& pts) {
    const int m = static_cast<int>(face.size());
    if (m < 3) return;
    if (m == 3) {
      triangles.push_back({pts[face[0]].index, pts[face[1]].index, pts[face[2]].index});
      return;
    }

    // Find top and bottom vertices by (y desc, x asc)
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

    // Find indices in face
    int idx_top = 0, idx_bot = 0;
    for (int i = 0; i < m; ++i) {
      if (face[i] == top) idx_top = i;
      if (face[i] == bottom) idx_bot = i;
    }

    // Build two chains from top to bottom
    std::vector<int> chain1, chain2;
    for (int i = idx_top; ; i = (i + 1) % m) {
      chain1.push_back(face[i]);
      if (face[i] == bottom) break;
    }
    for (int i = idx_top; ; i = (i - 1 + m) % m) {
      chain2.push_back(face[i]);
      if (face[i] == bottom) break;
    }

    // Determine left vs right by x-coordinate of second vertex
    double x_next1 = pts[chain1.size() > 1 ? chain1[1] : bottom].x;
    double x_next2 = pts[chain2.size() > 1 ? chain2[1] : bottom].x;
    bool chain1_is_left = x_next1 < x_next2;

    const auto& left_chain = chain1_is_left ? chain1 : chain2;
    const auto& right_chain = chain1_is_left ? chain2 : chain1;

    // Ensure buffer is large enough
    const int n = static_cast<int>(pts.size());
    if (static_cast<int>(is_left_buf_.size()) < n) {
      is_left_buf_.assign(n, 0);
    }
    // Mark left chain
    for (int v : is_left_touched_) is_left_buf_[v] = 0;
    is_left_touched_.clear();
    for (int v : left_chain) {
      if (!is_left_buf_[v]) {
        is_left_buf_[v] = 1;
        is_left_touched_.push_back(v);
      }
    }

    // Merge chains into sorted order (by y, top to bottom)
    std::vector<int> sorted;
    sorted.reserve(m);
    sorted.push_back(top);
    size_t i = 1, j = 1;
    while (i + 1 < left_chain.size() || j + 1 < right_chain.size()) {
      if (i + 1 >= left_chain.size()) {
        sorted.push_back(right_chain[j++]);
      } else if (j + 1 >= right_chain.size()) {
        sorted.push_back(left_chain[i++]);
      } else if (better_top(left_chain[i], right_chain[j])) {
        sorted.push_back(left_chain[i++]);
      } else {
        sorted.push_back(right_chain[j++]);
      }
    }
    sorted.push_back(bottom);

    // Stack-based monotone triangulation
    std::vector<int> st;
    st.push_back(sorted[0]);
    st.push_back(sorted[1]);

    for (int k = 2; k < m - 1; ++k) {
      int v = sorted[k];
      if (is_left_buf_[v] != is_left_buf_[st.back()]) {
        // Different chains
        while (static_cast<int>(st.size()) > 1) {
          int u = st.back(); st.pop_back();
          int w = st.back();
          triangles.push_back({pts[v].index, pts[u].index, pts[w].index});
        }
        st.pop_back();
        st.push_back(sorted[k - 1]);
        st.push_back(v);
      } else {
        // Same chain
        int u = st.back(); st.pop_back();
        while (!st.empty()) {
          int w = st.back();
          double c = detail::cross(pts[v], pts[u], pts[w]);
          bool ok = is_left_buf_[v] ? (c > detail::kEps) : (c < -detail::kEps);
          if (!ok) break;
          triangles.push_back({pts[v].index, pts[u].index, pts[w].index});
          u = st.back(); st.pop_back();
        }
        st.push_back(u);
        st.push_back(v);
      }
    }

    // Connect bottom to remaining stack
    int v = sorted[m - 1];
    while (static_cast<int>(st.size()) > 1) {
      int u = st.back(); st.pop_back();
      int w = st.back();
      triangles.push_back({pts[v].index, pts[u].index, pts[w].index});
    }
  }
};


}  // namespace chain_optimal
