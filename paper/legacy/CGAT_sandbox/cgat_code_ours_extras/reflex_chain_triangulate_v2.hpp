#pragma once
/**
 * Chain-based polygon triangulation: TRUE O(n + k log k) implementation.
 *
 * k = number of local maxima (equivalently, local minima) w.r.t. sweep direction.
 *
 * Complexity breakdown:
 * - Phase 1 (chain construction): O(n) - linear scan
 * - Phase 2 (monotone decomposition): O(k log k) - only 2k events, BST of size O(k)
 * - Phase 3 (triangulation): O(n) - half-edge face extraction + linear triangulation
 *
 * Key insight for O(n) face extraction:
 * Use a doubly-linked vertex list where adding a diagonal (u,v) duplicates
 * vertices and rewires pointers in O(1), splitting the polygon into two loops.
 * Each loop is then a monotone face that can be triangulated in O(face size).
 *
 * No heuristics. No fallbacks. Clean algorithm as described in the paper.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace reflex_chain_v2 {

struct Point {
  double x = 0.0;
  double y = 0.0;
  int index = -1;  // original vertex index
};

struct Triangle {
  int v0, v1, v2;
};

namespace detail {

static constexpr double kEps = 1e-12;

inline double cross(double ax, double ay, double bx, double by, double cx, double cy) {
  return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

inline double cross(const Point& a, const Point& b, const Point& c) {
  return cross(a.x, a.y, b.x, b.y, c.x, c.y);
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

// Strict total order: (y ascending, then x descending, then index descending)
// "below(a, b)" means a is strictly below b in this order.
inline bool below(const std::vector<Point>& pts, int a, int b) {
  if (pts[a].y < pts[b].y) return true;
  if (pts[a].y > pts[b].y) return false;
  if (pts[a].x > pts[b].x) return true;
  if (pts[a].x < pts[b].x) return false;
  return a > b;
}

inline bool above(const std::vector<Point>& pts, int a, int b) {
  return below(pts, b, a);
}

}  // namespace detail

class Triangulator {
public:
  std::vector<Triangle> triangles;
  int reflex_count = 0;
  int k_count = 0;  // number of local maxima (= number of local minima)

  void triangulate(std::vector<Point>& pts) {
    triangles.clear();
    reflex_count = 0;
    k_count = 0;

    const int n = static_cast<int>(pts.size());
    if (n < 3) return;

    // Ensure CCW orientation
    if (detail::signed_area(pts) < 0.0) {
      std::reverse(pts.begin(), pts.end());
      // Fix indices after reversal
      for (int i = 0; i < n; ++i) pts[i].index = i;
    }

    // Count reflex vertices and local maxima
    for (int i = 0; i < n; ++i) {
      const int p = (i - 1 + n) % n;
      const int nx = (i + 1) % n;
      const double c = detail::cross(pts[p], pts[i], pts[nx]);
      if (c < -detail::kEps) ++reflex_count;
      if (detail::below(pts, p, i) && detail::below(pts, nx, i)) ++k_count;
    }

    // Fast path: convex polygon -> fan triangulation O(n)
    if (reflex_count == 0) {
      triangles.reserve(n - 2);
      for (int i = 1; i < n - 1; ++i) {
        triangles.push_back({pts[0].index, pts[i].index, pts[i + 1].index});
      }
      return;
    }

    // Fast path: triangle
    if (n == 3) {
      triangles.push_back({pts[0].index, pts[1].index, pts[2].index});
      return;
    }

    // Main algorithm
    build_vertex_types(pts);
    build_chains(pts);
    init_half_edges(n);
    decompose_monotone(pts);
    extract_and_triangulate_faces(pts);
  }

private:
  enum class VType : uint8_t { Start, End, Split, Merge, Regular };

  std::vector<VType> types_;
  
  // Chain: sequence of vertices from local max to local min along CCW boundary
  struct Chain {
    std::vector<int> verts;  // vertex indices from upper to lower
    int curr = 0;            // current edge pointer for lazy advancement
    int pending = -1;        // pending merge vertex (helper)
    int upper = -1;          // local max vertex
    int lower = -1;          // local min vertex
  };
  std::vector<Chain> chains_;
  std::vector<int> max_to_chain_;  // local max -> left-boundary chain id
  std::vector<int> min_to_chain_;  // local min -> left-boundary chain id

  // Half-edge structure for O(1) diagonal insertion and O(n) face extraction
  struct HalfEdge {
    int target;    // vertex this half-edge points to
    int next;      // next half-edge in face (CCW)
    int twin;      // opposite half-edge
    int face;      // face id (-1 = unassigned)
  };
  std::vector<HalfEdge> hedges_;
  std::vector<int> v_hedge_;  // vertex -> one outgoing half-edge

  // Sweep state
  double sweep_y_ = 0.0;

  // Treap for status structure (chains ordered by x at sweep_y_)
  struct TreapNode {
    int chain_id;
    uint32_t prio;
    TreapNode* left = nullptr;
    TreapNode* right = nullptr;
  };
  TreapNode* status_root_ = nullptr;
  uint32_t rng_ = 0xDEADBEEF;

  uint32_t next_prio() {
    uint32_t x = rng_;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng_ = x;
    return x;
  }

  // ============ Phase 1: Classify vertices ============
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

      if (p_below && nx_below) {
        types_[i] = convex ? VType::Start : VType::Split;
      } else if (!p_below && !nx_below) {
        types_[i] = convex ? VType::End : VType::Merge;
      }
    }
  }

  // ============ Phase 1: Build monotone chains ============
  void build_chains(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());
    chains_.clear();
    max_to_chain_.assign(n, -1);
    min_to_chain_.assign(n, -1);

    for (int max_v = 0; max_v < n; ++max_v) {
      if (types_[max_v] != VType::Start && types_[max_v] != VType::Split) continue;

      Chain ch;
      ch.verts.push_back(max_v);
      ch.upper = max_v;

      int cur = (max_v + 1) % n;
      int safety = 0;
      while (types_[cur] != VType::End && types_[cur] != VType::Merge) {
        ch.verts.push_back(cur);
        cur = (cur + 1) % n;
        if (++safety > n + 2) throw std::runtime_error("chain loop");
      }
      ch.verts.push_back(cur);
      ch.lower = cur;
      ch.curr = 0;
      ch.pending = -1;

      int id = static_cast<int>(chains_.size());
      chains_.push_back(std::move(ch));
      max_to_chain_[max_v] = id;
      min_to_chain_[cur] = id;
    }
  }

  // ============ Half-edge initialization ============
  void init_half_edges(int n) {
    // Initialize half-edges for polygon boundary only
    // Diagonals will be added during sweep
    hedges_.clear();
    hedges_.reserve(4 * n);  // boundary + diagonals
    v_hedge_.assign(n, -1);

    // Create boundary half-edges
    for (int i = 0; i < n; ++i) {
      int j = (i + 1) % n;
      int he_ij = static_cast<int>(hedges_.size());
      hedges_.push_back({j, -1, -1, -1});
      int he_ji = static_cast<int>(hedges_.size());
      hedges_.push_back({i, -1, -1, -1});

      hedges_[he_ij].twin = he_ji;
      hedges_[he_ji].twin = he_ij;

      if (v_hedge_[i] < 0) v_hedge_[i] = he_ij;
    }

    // Set next pointers for boundary (initially one face)
    for (int i = 0; i < n; ++i) {
      int j = (i + 1) % n;
      int he_ij = 2 * i;      // half-edge i -> j
      int he_jk = 2 * j;      // half-edge j -> k where k = (j+1)%n
      hedges_[he_ij].next = he_jk;

      // Opposite direction: j -> i, next is prev boundary edge going backward
      int he_ji = 2 * i + 1;
      int prev_i = (i - 1 + n) % n;
      int he_prev = 2 * prev_i + 1;  // half-edge i -> prev_i
      hedges_[he_ji].next = he_prev;
    }
  }

  // Add diagonal (u, v) to half-edge structure in O(1)
  void add_diagonal(int u, int v, const std::vector<Point>& pts) {
    if (u == v) return;

    // Find half-edges leaving u and v that will be split by this diagonal
    // We need he_u such that the diagonal (u,v) lies between he_u.prev and he_u
    // This requires finding the correct angular position

    int he_uv = static_cast<int>(hedges_.size());
    hedges_.push_back({v, -1, -1, -1});
    int he_vu = static_cast<int>(hedges_.size());
    hedges_.push_back({u, -1, -1, -1});

    hedges_[he_uv].twin = he_vu;
    hedges_[he_vu].twin = he_uv;

    // Find the correct half-edges to rewire
    // For u: find outgoing half-edge he_u such that v is CCW from he_u's target
    // For v: find outgoing half-edge he_v such that u is CCW from he_v's target

    auto find_insert_pos = [&](int from, int to) -> int {
      int start = v_hedge_[from];
      if (start < 0) return -1;
      
      int best = start;
      int cur = start;
      do {
        int next_he = hedges_[hedges_[cur].twin].next;
        if (next_he < 0) break;
        
        int t1 = hedges_[cur].target;
        int t2 = hedges_[next_he].target;
        
        // Check if 'to' lies in the angular sector from t1 to t2 around 'from'
        double c1 = detail::cross(pts[from].x, pts[from].y,
                                   pts[t1].x, pts[t1].y,
                                   pts[to].x, pts[to].y);
        double c2 = detail::cross(pts[from].x, pts[from].y,
                                   pts[to].x, pts[to].y,
                                   pts[t2].x, pts[t2].y);
        double c12 = detail::cross(pts[from].x, pts[from].y,
                                    pts[t1].x, pts[t1].y,
                                    pts[t2].x, pts[t2].y);
        
        if (c12 >= 0) {
          // Convex sector: to must be left of both
          if (c1 >= 0 && c2 >= 0) {
            best = cur;
            break;
          }
        } else {
          // Reflex sector: to must be left of either
          if (c1 >= 0 || c2 >= 0) {
            best = cur;
            break;
          }
        }
        
        cur = next_he;
      } while (cur != start);
      
      return best;
    };

    int he_u = find_insert_pos(u, v);
    int he_v = find_insert_pos(v, u);

    if (he_u < 0 || he_v < 0) return;

    // Get the half-edges that currently follow he_u and he_v
    int he_u_next = hedges_[he_u].next;
    int he_v_next = hedges_[he_v].next;

    // Rewire: insert diagonal
    // Before: ... -> he_u -> he_u_next -> ...
    //         ... -> he_v -> he_v_next -> ...
    // After:  ... -> he_u -> he_uv -> he_v_next -> ...
    //         ... -> he_v -> he_vu -> he_u_next -> ...

    hedges_[he_u].next = he_uv;
    hedges_[he_uv].next = he_v_next;
    hedges_[he_v].next = he_vu;
    hedges_[he_vu].next = he_u_next;
  }

  // ============ Chain x-coordinate at sweep level ============
  void advance_chain(int cid, const std::vector<Point>& pts) {
    Chain& ch = chains_[cid];
    const int n = static_cast<int>(pts.size());

    while (ch.curr + 1 < static_cast<int>(ch.verts.size()) &&
           pts[ch.verts[ch.curr + 1]].y > sweep_y_ + detail::kEps) {
      ++ch.curr;
      int vreg = ch.verts[ch.curr];
      if (types_[vreg] != VType::Regular) continue;
      if (ch.pending < 0) continue;

      // Regular vertex with interior to right: connect to pending merge
      int prev = (vreg - 1 + n) % n;
      if (detail::below(pts, vreg, prev)) {
        add_diagonal(vreg, ch.pending, pts);
        ch.pending = -1;
      }
    }

    if (ch.curr + 1 >= static_cast<int>(ch.verts.size())) {
      ch.curr = static_cast<int>(ch.verts.size()) - 2;
    }
    if (ch.curr < 0) ch.curr = 0;
  }

  double chain_x_at(int cid, const std::vector<Point>& pts) {
    advance_chain(cid, pts);
    const Chain& ch = chains_[cid];
    int a = ch.verts[ch.curr];
    int b = ch.verts[ch.curr + 1];
    const auto& pa = pts[a];
    const auto& pb = pts[b];
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
    double xa = chain_x_at(a, pts);
    double xb = chain_x_at(b, pts);
    if (xa < xb) return true;
    if (xb < xa) return false;
    return a < b;
  }

  // ============ Treap operations ============
  TreapNode* treap_merge(TreapNode* a, TreapNode* b) {
    if (!a) return b;
    if (!b) return a;
    if (a->prio < b->prio) {
      a->right = treap_merge(a->right, b);
      return a;
    }
    b->left = treap_merge(a, b->left);
    return b;
  }

  void treap_split(TreapNode* root, int pivot, const std::vector<Point>& pts,
                   TreapNode*& left, TreapNode*& right) {
    if (!root) {
      left = right = nullptr;
      return;
    }
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
      if (chain_less(cid, (*cur)->chain_id, pts)) {
        cur = &(*cur)->left;
      } else {
        cur = &(*cur)->right;
      }
    }
  }

  int treap_predecessor_by_x(double x, const std::vector<Point>& pts) {
    int best = -1;
    TreapNode* cur = status_root_;
    while (cur) {
      double cx = chain_x_at(cur->chain_id, pts);
      if (cx < x) {
        best = cur->chain_id;
        cur = cur->right;
      } else {
        cur = cur->left;
      }
    }
    return best;
  }

  void treap_free(TreapNode* node) {
    if (!node) return;
    treap_free(node->left);
    treap_free(node->right);
    delete node;
  }

  // ============ Phase 2: Monotone decomposition O(k log k) ============
  void decompose_monotone(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());

    // Collect only extrema events: O(k) events
    struct Event {
      int v;
      double y, x;
      VType type;
    };
    std::vector<Event> events;
    events.reserve(2 * k_count);

    for (int i = 0; i < n; ++i) {
      if (types_[i] == VType::Start || types_[i] == VType::Split ||
          types_[i] == VType::End || types_[i] == VType::Merge) {
        events.push_back({i, pts[i].y, pts[i].x, types_[i]});
      }
    }

    // Sort events: O(k log k)
    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
      if (std::abs(a.y - b.y) > detail::kEps) return a.y > b.y;
      return a.x < b.x;
    });

    // Reset chain state
    for (auto& c : chains_) {
      c.curr = 0;
      c.pending = -1;
    }
    status_root_ = nullptr;

    // Process events: O(k) events, each O(log k)
    for (const auto& ev : events) {
      int v = ev.v;
      double vy = pts[v].y;

      // Set sweep_y_ appropriately
      if (ev.type == VType::End || ev.type == VType::Merge) {
        sweep_y_ = std::nextafter(vy, std::numeric_limits<double>::infinity());
      } else {
        sweep_y_ = std::nextafter(vy, -std::numeric_limits<double>::infinity());
      }

      switch (ev.type) {
        case VType::Start: {
          int cid = max_to_chain_[v];
          treap_insert(cid, pts);
          break;
        }

        case VType::End: {
          int cid = min_to_chain_[v];
          advance_chain(cid, pts);
          if (chains_[cid].pending >= 0) {
            add_diagonal(v, chains_[cid].pending, pts);
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
              add_diagonal(v, chains_[left_cid].pending, pts);
              chains_[left_cid].pending = -1;
            } else {
              int tgt = chain_slab_entry(left_cid, pts);
              add_diagonal(v, tgt, pts);
            }
          }
          int cid = max_to_chain_[v];
          treap_insert(cid, pts);
          break;
        }

        case VType::Merge: {
          int cid = min_to_chain_[v];
          advance_chain(cid, pts);
          if (chains_[cid].pending >= 0) {
            add_diagonal(v, chains_[cid].pending, pts);
            chains_[cid].pending = -1;
          }
          treap_erase(cid, pts);

          int left_cid = treap_predecessor_by_x(pts[v].x, pts);
          if (left_cid >= 0) {
            advance_chain(left_cid, pts);
            if (chains_[left_cid].pending >= 0) {
              add_diagonal(v, chains_[left_cid].pending, pts);
            }
            chains_[left_cid].pending = v;
          }
          break;
        }

        default:
          break;
      }
    }

    treap_free(status_root_);
    status_root_ = nullptr;
  }

  // ============ Phase 3: Face extraction + triangulation O(n) ============
  void extract_and_triangulate_faces(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());
    const int m = static_cast<int>(hedges_.size());

    // Mark all half-edges as unvisited
    std::vector<bool> visited(m, false);

    // Extract faces by walking half-edge cycles
    for (int start = 0; start < m; ++start) {
      if (visited[start]) continue;

      // Walk the face
      std::vector<int> face_verts;
      int he = start;
      int safety = 0;
      do {
        if (he < 0 || he >= m || safety > m) break;
        visited[he] = true;
        // The vertex at the START of this half-edge
        int twin = hedges_[he].twin;
        if (twin >= 0 && twin < m) {
          face_verts.push_back(hedges_[twin].target);
        }
        he = hedges_[he].next;
        ++safety;
      } while (he != start && he >= 0);

      if (face_verts.size() < 3) continue;

      // Compute signed area to identify interior vs exterior face
      double area = 0.0;
      for (size_t i = 0; i < face_verts.size(); ++i) {
        size_t j = (i + 1) % face_verts.size();
        int vi = face_verts[i];
        int vj = face_verts[j];
        area += pts[vi].x * pts[vj].y - pts[vj].x * pts[vi].y;
      }

      // Interior faces have positive area (CCW), exterior has negative (CW)
      if (area > detail::kEps) {
        triangulate_monotone_face(face_verts, pts);
      }
    }
  }

  // Triangulate a y-monotone face in O(face size)
  void triangulate_monotone_face(const std::vector<int>& face, const std::vector<Point>& pts) {
    const int m = static_cast<int>(face.size());
    if (m < 3) return;
    if (m == 3) {
      triangles.push_back({pts[face[0]].index, pts[face[1]].index, pts[face[2]].index});
      return;
    }

    // Find topmost and bottommost vertices
    int top_idx = 0, bot_idx = 0;
    for (int i = 1; i < m; ++i) {
      if (detail::above(pts, face[i], face[top_idx])) top_idx = i;
      if (detail::below(pts, face[i], face[bot_idx])) bot_idx = i;
    }

    // Build left and right chains from top to bottom
    std::vector<int> left_chain, right_chain;
    
    // Walk CCW from top_idx to bot_idx -> this is the "left" chain
    int i = top_idx;
    while (i != bot_idx) {
      left_chain.push_back(face[i]);
      i = (i + 1) % m;
    }
    left_chain.push_back(face[bot_idx]);

    // Walk CCW from bot_idx to top_idx -> this is the "right" chain (reversed)
    i = bot_idx;
    while (i != top_idx) {
      right_chain.push_back(face[i]);
      i = (i + 1) % m;
    }
    right_chain.push_back(face[top_idx]);
    std::reverse(right_chain.begin(), right_chain.end());

    // Merge chains into sorted vertex list with chain membership
    std::vector<std::pair<int, bool>> sorted;  // (vertex, is_left_chain)
    sorted.reserve(m);
    
    size_t li = 0, ri = 0;
    while (li < left_chain.size() || ri < right_chain.size()) {
      if (li >= left_chain.size()) {
        sorted.push_back({right_chain[ri++], false});
      } else if (ri >= right_chain.size()) {
        sorted.push_back({left_chain[li++], true});
      } else if (detail::above(pts, left_chain[li], right_chain[ri])) {
        sorted.push_back({left_chain[li++], true});
      } else {
        sorted.push_back({right_chain[ri++], false});
      }
    }

    // Standard monotone polygon triangulation
    std::vector<std::pair<int, bool>> stack;
    stack.push_back(sorted[0]);
    stack.push_back(sorted[1]);

    for (size_t j = 2; j < sorted.size(); ++j) {
      int vj = sorted[j].first;
      bool side_j = sorted[j].second;

      if (side_j != stack.back().second) {
        // Different chain: triangulate with all stack vertices
        while (stack.size() > 1) {
          int v1 = stack[stack.size() - 1].first;
          int v2 = stack[stack.size() - 2].first;
          triangles.push_back({pts[vj].index, pts[v1].index, pts[v2].index});
          stack.pop_back();
        }
        stack.pop_back();
        stack.push_back(sorted[j - 1]);
        stack.push_back(sorted[j]);
      } else {
        // Same chain: pop while diagonal is inside
        auto last = stack.back();
        stack.pop_back();
        while (!stack.empty()) {
          int v1 = last.first;
          int v2 = stack.back().first;
          double c = detail::cross(pts[vj], pts[v1], pts[v2]);
          bool inside = side_j ? (c > detail::kEps) : (c < -detail::kEps);
          if (!inside) break;
          triangles.push_back({pts[vj].index, pts[v1].index, pts[v2].index});
          last = stack.back();
          stack.pop_back();
        }
        stack.push_back(last);
        stack.push_back(sorted[j]);
      }
    }
  }
};

}  // namespace reflex_chain_v2
