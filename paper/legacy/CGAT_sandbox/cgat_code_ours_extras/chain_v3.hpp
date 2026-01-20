#pragma once
/**
 * Chain-based polygon triangulation: O(n + k log k)
 *
 * This implementation:
 * 1. Uses chain sweep for O(k log k) diagonal generation
 * 2. Collects diagonals into a list
 * 3. Uses PolyPartition-style algorithm for face extraction (O(n log n) but fast in practice)
 *
 * Phase 1: O(n) - vertex classification, chain construction
 * Phase 2: O(k log k) - sweep with k events, BST of O(k) chains
 * Phase 3: O(n log n) - face extraction using angle sorting (bottleneck)
 *
 * True O(n + k log k) requires half-edge structure with O(1) diagonal insertion,
 * which is complex to implement correctly. This version is a practical compromise.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include <vector>

namespace chain_v3 {

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

    if (detail::signed_area(pts) < 0.0) {
      std::reverse(pts.begin(), pts.end());
      for (int i = 0; i < n; ++i) pts[i].index = i;
    }

    for (int i = 0; i < n; ++i) {
      const int p = (i - 1 + n) % n;
      const int nx = (i + 1) % n;
      if (detail::cross(pts[p], pts[i], pts[nx]) < -detail::kEps) ++reflex_count;
      if (detail::below(pts, p, i) && detail::below(pts, nx, i)) ++k_count;
    }

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

    build_vertex_types(pts);
    build_chains(pts);
    generate_diagonals(pts);
    
    diagonal_count = static_cast<int>(diagonals_.size());
    
    extract_and_triangulate_faces(pts);
  }

private:
  enum class VType : uint8_t { Start, End, Split, Merge, Regular };
  
  std::vector<VType> types_;
  std::vector<std::pair<int, int>> diagonals_;
  std::unordered_set<uint64_t> diag_set_;

  struct Chain {
    std::vector<int> verts;
    int curr = 0;
    int pending = -1;
  };
  std::vector<Chain> chains_;
  std::vector<int> max_to_chain_;
  std::vector<int> min_to_chain_;

  struct TreapNode {
    int chain_id;
    uint32_t prio;
    TreapNode* left = nullptr;
    TreapNode* right = nullptr;
  };
  TreapNode* status_root_ = nullptr;
  uint32_t rng_ = 0xFEEDBEEF;
  double sweep_y_ = 0.0;

  uint32_t next_prio() {
    uint32_t x = rng_;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng_ = x;
    return x;
  }

  static uint64_t diag_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
  }

  void add_diagonal(int u, int v) {
    if (u == v) return;
    uint64_t key = diag_key(u, v);
    if (diag_set_.insert(key).second) {
      diagonals_.push_back({std::min(u, v), std::max(u, v)});
    }
  }

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

  // Phase 3: face extraction + triangulation (using adjacency sorting)
  void extract_and_triangulate_faces(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());

    // Build adjacency
    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; ++i) {
      adj[i].push_back((i + 1) % n);
      adj[(i + 1) % n].push_back(i);
    }
    for (const auto& [a, b] : diagonals_) {
      adj[a].push_back(b);
      adj[b].push_back(a);
    }

    // Sort neighbors by angle (CCW)
    for (int u = 0; u < n; ++u) {
      auto& nb = adj[u];
      std::sort(nb.begin(), nb.end());
      nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
      if (nb.size() <= 2) continue;
      std::sort(nb.begin(), nb.end(), [&](int a, int b) {
        double ax = pts[a].x - pts[u].x, ay = pts[a].y - pts[u].y;
        double bx = pts[b].x - pts[u].x, by = pts[b].y - pts[u].y;
        int ha = (ay > 0 || (ay == 0 && ax >= 0)) ? 0 : 1;
        int hb = (by > 0 || (by == 0 && bx >= 0)) ? 0 : 1;
        if (ha != hb) return ha < hb;
        double cr = ax * by - ay * bx;
        if (cr > 0) return true;
        if (cr < 0) return false;
        return ax * ax + ay * ay < bx * bx + by * by;
      });
    }

    // Build offset array for half-edge indexing
    std::vector<int> off(n + 1, 0);
    for (int u = 0; u < n; ++u) off[u + 1] = off[u] + static_cast<int>(adj[u].size());
    const int m_dir = off[n];

    // Reverse position lookup
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
          if (prev_idx < 0) break;
          int next_i = (prev_idx - 1 + static_cast<int>(adj[curr].size())) % static_cast<int>(adj[curr].size());
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

    // Identify outer face (largest absolute area) and triangulate interior faces
    if (faces.empty()) return;
    size_t outer_idx = 0;
    double best_abs = std::abs(face_areas[0]);
    for (size_t i = 1; i < faces.size(); ++i) {
      if (std::abs(face_areas[i]) > best_abs) {
        best_abs = std::abs(face_areas[i]);
        outer_idx = i;
      }
    }

    for (size_t i = 0; i < faces.size(); ++i) {
      if (i == outer_idx) continue;
      auto& face = faces[i];
      if (face_areas[i] < 0) std::reverse(face.begin(), face.end());
      triangulate_monotone_face(face, pts);
    }
  }

  void triangulate_monotone_face(const std::vector<int>& face, const std::vector<Point>& pts) {
    const int m = static_cast<int>(face.size());
    if (m < 3) return;
    if (m == 3) {
      triangles.push_back({pts[face[0]].index, pts[face[1]].index, pts[face[2]].index});
      return;
    }

    // Find top/bottom
    auto above = [&](int a, int b) {
      if (pts[a].y > pts[b].y) return true;
      if (pts[a].y < pts[b].y) return false;
      return pts[a].x < pts[b].x;
    };
    int top = 0, bot = 0;
    for (int i = 1; i < m; ++i) {
      if (above(face[i], face[top])) top = i;
      if (!above(face[i], face[bot])) bot = i;
    }

    // Build chains
    std::vector<int> left_chain, right_chain;
    for (int i = top; ; i = (i + 1) % m) {
      left_chain.push_back(i);
      if (i == bot) break;
    }
    for (int i = bot; ; i = (i + 1) % m) {
      right_chain.push_back(i);
      if (i == top) break;
    }
    std::reverse(right_chain.begin(), right_chain.end());

    // Merge sorted
    std::vector<std::pair<int, bool>> sorted;
    size_t li = 0, ri = 0;
    while (li < left_chain.size() || ri < right_chain.size()) {
      if (li >= left_chain.size()) sorted.push_back({right_chain[ri++], false});
      else if (ri >= right_chain.size()) sorted.push_back({left_chain[li++], true});
      else if (above(face[left_chain[li]], face[right_chain[ri]])) sorted.push_back({left_chain[li++], true});
      else sorted.push_back({right_chain[ri++], false});
    }

    // Standard monotone triangulation
    std::vector<std::pair<int, bool>> stack;
    stack.push_back(sorted[0]);
    stack.push_back(sorted[1]);

    for (size_t j = 2; j < sorted.size(); ++j) {
      int vj = sorted[j].first;
      bool side_j = sorted[j].second;
      if (side_j != stack.back().second) {
        while (stack.size() > 1) {
          int v1 = stack.back().first; stack.pop_back();
          int v2 = stack.back().first;
          triangles.push_back({pts[face[vj]].index, pts[face[v1]].index, pts[face[v2]].index});
        }
        stack.clear();
        stack.push_back(sorted[j - 1]);
        stack.push_back(sorted[j]);
      } else {
        auto last = stack.back(); stack.pop_back();
        while (!stack.empty()) {
          int v1 = last.first, v2 = stack.back().first;
          double c = detail::cross(pts[face[vj]], pts[face[v1]], pts[face[v2]]);
          bool inside = side_j ? (c > detail::kEps) : (c < -detail::kEps);
          if (!inside) break;
          triangles.push_back({pts[face[vj]].index, pts[face[v1]].index, pts[face[v2]].index});
          last = stack.back(); stack.pop_back();
        }
        stack.push_back(last);
        stack.push_back(sorted[j]);
      }
    }
  }
};

}  // namespace chain_v3
