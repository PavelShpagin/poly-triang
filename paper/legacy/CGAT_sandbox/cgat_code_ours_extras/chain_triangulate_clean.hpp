#pragma once
/**
 * Chain-based polygon triangulation: O(n + k log k)
 *
 * k = number of local maxima (= local minima) w.r.t. sweep direction.
 *
 * Phase 1 (O(n)): Classify vertices, build chains
 * Phase 2 (O(k log k)): Sweep to generate diagonals (2k events, BST of O(k) chains)
 * Phase 3 (O(n)): PolyPartition-style linked-list face splitting + monotone triangulation
 *
 * NO heuristics. NO fallbacks. Clean implementation.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace chain_clean {

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

inline bool above(const std::vector<Point>& pts, int a, int b) {
  return below(pts, b, a);
}
}  // namespace detail

class Triangulator {
public:
  std::vector<Triangle> triangles;
  int reflex_count = 0;
  int k_count = 0;
  int diagonal_count = 0;
  int debug_face_count = 0;
  int debug_interior_faces = 0;

  void triangulate(std::vector<Point>& pts) {
    triangles.clear();
    reflex_count = 0;
    k_count = 0;
    diagonal_count = 0;

    const int n = static_cast<int>(pts.size());
    if (n < 3) return;

    // Ensure CCW
    if (detail::signed_area(pts) < 0.0) {
      std::reverse(pts.begin(), pts.end());
      for (int i = 0; i < n; ++i) pts[i].index = i;
    }

    // Count reflex and k
    for (int i = 0; i < n; ++i) {
      const int p = (i - 1 + n) % n;
      const int nx = (i + 1) % n;
      if (detail::cross(pts[p], pts[i], pts[nx]) < -detail::kEps) ++reflex_count;
      if (detail::below(pts, p, i) && detail::below(pts, nx, i)) ++k_count;
    }

    // Fast path: convex
    if (reflex_count == 0) {
      triangles.reserve(n - 2);
      for (int i = 1; i < n - 1; ++i) {
        triangles.push_back({pts[0].index, pts[i].index, pts[i + 1].index});
      }
      return;
    }

    if (n == 3) {
      triangles.push_back({pts[0].index, pts[1].index, pts[2].index});
      return;
    }

    // Initialize linked list
    init_linked_list(pts);
    
    // Build chains and vertex types
    build_vertex_types(pts);
    build_chains(pts);
    
    // Generate diagonals via chain sweep
    generate_diagonals(pts);
    
    // Extract and triangulate faces
    triangulate_all_faces(pts);
  }

private:
  enum class VType : uint8_t { Start, End, Split, Merge, Regular };

  // Linked list vertex
  struct Vertex {
    double x, y;
    int orig_idx;
    int prev, next;
  };
  std::vector<Vertex> verts_;
  int num_verts_;

  std::vector<VType> types_;

  // Chain
  struct Chain {
    std::vector<int> verts;
    int curr = 0;
    int pending = -1;
    int upper = -1;
    int lower = -1;
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
  uint32_t rng_ = 0xDEADC0DE;
  double sweep_y_ = 0.0;

  uint32_t next_prio() {
    uint32_t x = rng_;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng_ = x;
    return x;
  }

  void init_linked_list(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());
    verts_.resize(3 * n);
    num_verts_ = n;
    for (int i = 0; i < n; ++i) {
      verts_[i].x = pts[i].x;
      verts_[i].y = pts[i].y;
      verts_[i].orig_idx = pts[i].index;
      verts_[i].prev = (i - 1 + n) % n;
      verts_[i].next = (i + 1) % n;
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
      if (p_below && nx_below) {
        types_[i] = convex ? VType::Start : VType::Split;
      } else if (!p_below && !nx_below) {
        types_[i] = convex ? VType::End : VType::Merge;
      }
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
      int id = static_cast<int>(chains_.size());
      chains_.push_back(std::move(ch));
      max_to_chain_[max_v] = id;
      min_to_chain_[cur] = id;
    }
  }

  // Add diagonal using PolyPartition-style linked list split
  void add_diagonal(int index1, int index2) {
    if (index1 == index2) return;
    if (index1 >= num_verts_ || index2 >= num_verts_) return;

    int newindex1 = num_verts_++;
    int newindex2 = num_verts_++;

    verts_[newindex1] = verts_[index1];
    verts_[newindex2] = verts_[index2];

    verts_[newindex2].next = verts_[index2].next;
    verts_[newindex1].next = verts_[index1].next;

    verts_[verts_[index2].next].prev = newindex2;
    verts_[verts_[index1].next].prev = newindex1;

    verts_[index1].next = newindex2;
    verts_[newindex2].prev = index1;

    verts_[index2].next = newindex1;
    verts_[newindex1].prev = index2;

    ++diagonal_count;
  }

  void advance_chain(int cid, const std::vector<Point>& pts) {
    Chain& ch = chains_[cid];
    const int n = static_cast<int>(pts.size());
    while (ch.curr + 1 < static_cast<int>(ch.verts.size()) &&
           pts[ch.verts[ch.curr + 1]].y > sweep_y_ + detail::kEps) {
      ++ch.curr;
      int vreg = ch.verts[ch.curr];
      if (types_[vreg] != VType::Regular) continue;
      if (ch.pending < 0) continue;
      int prev = (vreg - 1 + n) % n;
      if (detail::below(pts, vreg, prev)) {
        add_diagonal(vreg, ch.pending);
        ch.pending = -1;
      }
    }
    if (ch.curr + 1 >= static_cast<int>(ch.verts.size())) {
      ch.curr = std::max(0, static_cast<int>(ch.verts.size()) - 2);
    }
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
      if (chain_less(cid, (*cur)->chain_id, pts)) cur = &(*cur)->left;
      else cur = &(*cur)->right;
    }
  }

  int treap_predecessor_by_x(double x, const std::vector<Point>& pts) {
    int best = -1;
    TreapNode* cur = status_root_;
    while (cur) {
      double cx = chain_x_at(cur->chain_id, pts);
      if (cx < x) { best = cur->chain_id; cur = cur->right; }
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
      if (types_[i] == VType::Start || types_[i] == VType::Split ||
          types_[i] == VType::End || types_[i] == VType::Merge) {
        events.push_back({i, pts[i].y, pts[i].x, types_[i]});
      }
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
      if (ev.type == VType::End || ev.type == VType::Merge)
        sweep_y_ = std::nextafter(vy, std::numeric_limits<double>::infinity());
      else
        sweep_y_ = std::nextafter(vy, -std::numeric_limits<double>::infinity());

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
            if (chains_[left_cid].pending >= 0) {
              add_diagonal(v, chains_[left_cid].pending);
            }
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

  void triangulate_all_faces(const std::vector<Point>& pts) {
    std::vector<bool> visited(num_verts_, false);
    debug_face_count = 0;
    debug_interior_faces = 0;
    
    for (int start = 0; start < num_verts_; ++start) {
      if (visited[start]) continue;
      std::vector<int> face;
      int cur = start;
      int safety = 0;
      while (!visited[cur] && safety < num_verts_ + 10) {
        visited[cur] = true;
        face.push_back(cur);
        cur = verts_[cur].next;
        ++safety;
      }
      if (face.size() < 3) continue;
      
      ++debug_face_count;
      
      // Compute signed area
      double area = 0.0;
      for (size_t i = 0; i < face.size(); ++i) {
        size_t j = (i + 1) % face.size();
        area += verts_[face[i]].x * verts_[face[j]].y -
                verts_[face[j]].x * verts_[face[i]].y;
      }
      // Interior faces have positive area (CCW)
      if (area > detail::kEps) {
        ++debug_interior_faces;
        triangulate_monotone_face(face, pts);
      }
    }
  }

  void triangulate_monotone_face(const std::vector<int>& face, const std::vector<Point>& pts) {
    const int m = static_cast<int>(face.size());
    if (m < 3) return;
    if (m == 3) {
      triangles.push_back({verts_[face[0]].orig_idx, verts_[face[1]].orig_idx, verts_[face[2]].orig_idx});
      return;
    }

    // Build coordinate array
    struct FV { double x, y; int orig; };
    std::vector<FV> fv(m);
    for (int i = 0; i < m; ++i) {
      fv[i].x = verts_[face[i]].x;
      fv[i].y = verts_[face[i]].y;
      fv[i].orig = verts_[face[i]].orig_idx;
    }

    // Find top/bottom
    int top = 0, bot = 0;
    for (int i = 1; i < m; ++i) {
      if (fv[i].y > fv[top].y || (fv[i].y == fv[top].y && fv[i].x < fv[top].x)) top = i;
      if (fv[i].y < fv[bot].y || (fv[i].y == fv[bot].y && fv[i].x > fv[bot].x)) bot = i;
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

    // Merge
    auto above_fv = [&](int a, int b) {
      if (fv[a].y > fv[b].y) return true;
      if (fv[a].y < fv[b].y) return false;
      return fv[a].x < fv[b].x;
    };
    std::vector<std::pair<int, bool>> sorted;
    size_t li = 0, ri = 0;
    while (li < left_chain.size() || ri < right_chain.size()) {
      if (li >= left_chain.size()) sorted.push_back({right_chain[ri++], false});
      else if (ri >= right_chain.size()) sorted.push_back({left_chain[li++], true});
      else if (above_fv(left_chain[li], right_chain[ri])) sorted.push_back({left_chain[li++], true});
      else sorted.push_back({right_chain[ri++], false});
    }

    // Standard monotone triangulation
    std::vector<std::pair<int, bool>> stack;
    stack.push_back(sorted[0]);
    stack.push_back(sorted[1]);

    auto cross_fv = [&](int a, int b, int c) {
      return (fv[b].x - fv[a].x) * (fv[c].y - fv[a].y) -
             (fv[b].y - fv[a].y) * (fv[c].x - fv[a].x);
    };

    for (size_t j = 2; j < sorted.size(); ++j) {
      int vj = sorted[j].first;
      bool side_j = sorted[j].second;
      if (side_j != stack.back().second) {
        while (stack.size() > 1) {
          int v1 = stack.back().first; stack.pop_back();
          int v2 = stack.back().first;
          triangles.push_back({fv[vj].orig, fv[v1].orig, fv[v2].orig});
        }
        stack.clear();
        stack.push_back(sorted[j - 1]);
        stack.push_back(sorted[j]);
      } else {
        auto last = stack.back(); stack.pop_back();
        while (!stack.empty()) {
          int v1 = last.first;
          int v2 = stack.back().first;
          double c = cross_fv(vj, v1, v2);
          bool inside = side_j ? (c > detail::kEps) : (c < -detail::kEps);
          if (!inside) break;
          triangles.push_back({fv[vj].orig, fv[v1].orig, fv[v2].orig});
          last = stack.back(); stack.pop_back();
        }
        stack.push_back(last);
        stack.push_back(sorted[j]);
      }
    }
  }
};

}  // namespace chain_clean
