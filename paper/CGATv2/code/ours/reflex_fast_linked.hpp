#pragma once
/**
 * Edge-based sweep-line triangulation with linked-list face splitting.
 *
 * This follows the classic y-monotone decomposition + triangulation pipeline
 * (PolyPartition-style). Inserted diagonals duplicate vertices and rewire a
 * doubly-linked list so polygon loops split in O(1), making face extraction fast.
 *
 * Used as a practical fallback in `reflex_cli.cpp` for high-k inputs.
 *
 * Worst-case time is O(n log n) due to sweep-line balanced-tree operations.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

namespace fast_linked {

struct Point {
  double x = 0.0;
  double y = 0.0;
  int index = -1;
};

struct Triangle {
  int v0, v1, v2;
};

class Triangulator {
public:
  std::vector<Triangle> triangles;
  int reflex_count = 0;
  int num_faces = 0;  // Debug: count faces
  int num_diagonals = 0;  // Debug: count diagonals

  void triangulate(const std::vector<Point>& input) {
    triangles.clear();
    reflex_count = 0;
    
    n_ = static_cast<int>(input.size());
    if (n_ < 3) return;
    
    pts_ = &input;
    
    // Check orientation
    double area = 0;
    for (int i = 0; i < n_; ++i) {
      int j = (i + 1) % n_;
      area += input[i].x * input[j].y - input[j].x * input[i].y;
    }
    if (area < 0) {
      reversed_ = true;
    } else {
      reversed_ = false;
    }
    
    // Count reflex vertices
    for (int i = 0; i < n_; ++i) {
      if (is_reflex(i)) ++reflex_count;
    }
    
    // Fast path for convex polygons
    if (reflex_count == 0) {
      triangles.reserve(static_cast<std::size_t>(n_ - 2));
      for (int i = 1; i < n_ - 1; ++i) {
        triangles.push_back({0, i, i + 1});
      }
      return;
    }
    
    if (n_ == 3) {
      triangles.push_back({0, 1, 2});
      return;
    }
    
    decompose_and_triangulate(input);
  }

private:
  static constexpr double kEps = 1e-12;
  
  int n_ = 0;
  const std::vector<Point>* pts_ = nullptr;
  bool reversed_ = false;
  
  // Linked-list vertex structure
  struct Vertex {
    double x, y;
    int next, prev;
  };
  std::vector<Vertex> verts_;
  int num_verts_ = 0;
  
  enum VType : uint8_t { START, END, SPLIT, MERGE, REGULAR };
  std::vector<VType> types_;
  
  // Edge in the sweep-line tree (matches PolyPartition's ScanLineEdge)
  struct Edge {
    mutable int index;
    mutable int helper;
    
    double x1, y1, x2, y2;  // Edge endpoints (p1 and p2)
    
    // IsConvex: returns true if p3 is to the left of line p1->p2
    static bool is_convex(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y) {
      double tmp = (p3y - p1y) * (p2x - p1x) - (p3x - p1x) * (p2y - p1y);
      return tmp > 0;
    }
    
    bool operator<(const Edge& o) const {
      // PolyPartition's comparison logic
      if (o.y1 == o.y2) {
        if (y1 == y2) {
          return y1 < o.y1;
        }
        return is_convex(x1, y1, x2, y2, o.x1, o.y1);
      } else if (y1 == y2) {
        return !is_convex(o.x1, o.y1, o.x2, o.y2, x1, y1);
      } else if (y1 < o.y1) {
        return !is_convex(o.x1, o.y1, o.x2, o.y2, x1, y1);
      } else {
        return is_convex(x1, y1, x2, y2, o.x1, o.y1);
      }
    }
  };
  
  std::set<Edge> edge_tree_;
  std::vector<std::set<Edge>::iterator> edge_iters_;
  
  int nxt(int i) const {
    return reversed_ ? ((i - 1 + n_) % n_) : ((i + 1) % n_);
  }
  
  int prv(int i) const {
    return reversed_ ? ((i + 1) % n_) : ((i - 1 + n_) % n_);
  }
  
  const Point& pt(int i) const {
    return (*pts_)[i];
  }
  
  double cross(int o, int a, int b) const {
    return (pt(a).x - pt(o).x) * (pt(b).y - pt(o).y) -
           (pt(a).y - pt(o).y) * (pt(b).x - pt(o).x);
  }
  
  bool is_reflex(int i) const {
    int p = prv(i), nx = nxt(i);
    double c = cross(p, i, nx);
    return c < -kEps;
  }
  
  bool below(int a, int b) const {
    if (std::abs(pt(a).y - pt(b).y) > kEps) return pt(a).y < pt(b).y;
    return pt(a).x > pt(b).x;
  }
  
  // PolyPartition-style diagonal insertion with linked-list rewiring
  // (From context summary)
  void add_diagonal(int i1, int i2) {
    if (i1 == i2) return;
    num_diagonals++;
    
    int new1 = num_verts_++;
    int new2 = num_verts_++;

    // Copy vertex data
    verts_[new1] = verts_[i1];
    verts_[new2] = verts_[i2];

    // Rewire the linked list - this splits it into two loops!
    verts_[new2].next = verts_[i2].next;
    verts_[new1].next = verts_[i1].next;

    verts_[verts_[i2].next].prev = new2;
    verts_[verts_[i1].next].prev = new1;

    verts_[i1].next = new2;
    verts_[new2].prev = i1;

    verts_[i2].next = new1;
    verts_[new1].prev = i2;

    // Copy types
    types_[new1] = types_[i1];
    types_[new2] = types_[i2];

    // Transfer edge iterators - new vertices take over the tree slots
    edge_iters_[new1] = edge_iters_[i1];
    edge_iters_[new2] = edge_iters_[i2];
    
    // Update the edge's stored index and helper to point to the new vertex
    if (edge_iters_[new1] != edge_tree_.end()) {
      edge_iters_[new1]->index = new1;
      edge_iters_[new1]->helper = edge_iters_[i1]->helper;
    }
    if (edge_iters_[new2] != edge_tree_.end()) {
      edge_iters_[new2]->index = new2;
      edge_iters_[new2]->helper = edge_iters_[i2]->helper;
    }
  }
  
  void decompose_and_triangulate(const std::vector<Point>& pts) {
    int maxverts = 3 * n_;
    verts_.resize(maxverts);
    types_.resize(maxverts);
    edge_iters_.resize(maxverts);
    num_verts_ = n_;
    
    // Initialize linked list
    for (int i = 0; i < n_; ++i) {
      verts_[i].x = pt(i).x;
      verts_[i].y = pt(i).y;
      verts_[i].next = nxt(i);
      verts_[i].prev = prv(i);
    }
    
    // Classify vertices
    for (int i = 0; i < n_; ++i) {
      int p = prv(i), nx = nxt(i);
      bool p_below = below(p, i);
      bool n_below = below(nx, i);
      double c = cross(p, i, nx);
      
      if (p_below && n_below) {
        types_[i] = (c > kEps) ? START : SPLIT;
      } else if (!p_below && !n_below) {
        types_[i] = (c > kEps) ? END : MERGE;
      } else {
        types_[i] = REGULAR;
      }
    }
    
    // Sort vertices by decreasing y
    std::vector<int> sorted(n_);
    for (int i = 0; i < n_; ++i) sorted[i] = i;
    std::sort(sorted.begin(), sorted.end(), [this](int a, int b) {
      if (std::abs(pt(a).y - pt(b).y) > kEps) return pt(a).y > pt(b).y;
      return pt(a).x < pt(b).x;
    });
    
    // Initialize edge iterators
    edge_tree_.clear();
    for (int i = 0; i < maxverts; ++i) {
      edge_iters_[i] = edge_tree_.end();
    }
    
    // Helpers array
    std::vector<int> helpers(maxverts, -1);
    
    // Process vertices in sweep order
    for (int vi : sorted) {
      int v = vi;
      int v2 = vi;
      
      switch (types_[vi]) {
        case START: {
          // Insert edge starting at v
          int vnext = verts_[v].next;
          Edge e{v, v, verts_[v].x, verts_[v].y, verts_[vnext].x, verts_[vnext].y};
          auto [it, ok] = edge_tree_.insert(e);
          edge_iters_[v] = it;
          helpers[v] = v;
          break;
        }
        
        case END: {
          int ep = verts_[v].prev;
          if (edge_iters_[ep] != edge_tree_.end()) {
            if (types_[helpers[ep]] == MERGE) {
              add_diagonal(v, helpers[ep]);
            }
            edge_tree_.erase(edge_iters_[ep]);
            edge_iters_[ep] = edge_tree_.end();
          }
          break;
        }
        
        case SPLIT: {
          // Find edge directly left
          Edge probe{v, v, verts_[v].x, verts_[v].y, verts_[v].x, verts_[v].y};
          auto it = edge_tree_.lower_bound(probe);
          if (it != edge_tree_.begin()) {
            --it;
            add_diagonal(v, helpers[it->index]);
            v2 = num_verts_ - 2;
            helpers[it->index] = v;
          }
          // Insert new edge
          int vnext = verts_[v2].next;
          Edge e{v2, v2, verts_[v2].x, verts_[v2].y, verts_[vnext].x, verts_[vnext].y};
          auto [eit, ok] = edge_tree_.insert(e);
          edge_iters_[v2] = eit;
          helpers[v2] = v2;
          break;
        }
        
        case MERGE: {
          int ep = verts_[v].prev;
          if (edge_iters_[ep] != edge_tree_.end()) {
            if (types_[helpers[ep]] == MERGE) {
              add_diagonal(v, helpers[ep]);
              v2 = num_verts_ - 2;
            }
            edge_tree_.erase(edge_iters_[ep]);
            edge_iters_[ep] = edge_tree_.end();
          }
          // Find edge directly left
          Edge probe{v2, v2, verts_[v2].x, verts_[v2].y, verts_[v2].x, verts_[v2].y};
          auto it = edge_tree_.lower_bound(probe);
          if (it != edge_tree_.begin()) {
            --it;
            if (types_[helpers[it->index]] == MERGE) {
              add_diagonal(v2, helpers[it->index]);
            }
            helpers[it->index] = v2;
          }
          break;
        }
        
        case REGULAR: {
          int ep = verts_[v].prev;
          // Check if interior is to the right
          if (verts_[ep].y > verts_[v].y + kEps ||
              (std::abs(verts_[ep].y - verts_[v].y) < kEps && verts_[ep].x < verts_[v].x)) {
            // Interior to the right
            if (edge_iters_[ep] != edge_tree_.end()) {
              if (types_[helpers[ep]] == MERGE) {
                add_diagonal(v, helpers[ep]);
                v2 = num_verts_ - 2;
              }
              edge_tree_.erase(edge_iters_[ep]);
              edge_iters_[ep] = edge_tree_.end();
            }
            int vnext = verts_[v2].next;
            Edge e{v2, v, verts_[v2].x, verts_[v2].y, verts_[vnext].x, verts_[vnext].y};
            auto [eit, ok] = edge_tree_.insert(e);
            edge_iters_[v2] = eit;
            helpers[v2] = v;
          } else {
            // Interior to the left
            Edge probe{v, v, verts_[v].x, verts_[v].y, verts_[v].x, verts_[v].y};
            auto it = edge_tree_.lower_bound(probe);
            if (it != edge_tree_.begin()) {
              --it;
              if (types_[helpers[it->index]] == MERGE) {
                add_diagonal(v, helpers[it->index]);
              }
              helpers[it->index] = v;
            }
          }
          break;
        }
      }
    }
    
    // Extract and triangulate monotone faces
    // Use character array like PolyPartition does
    std::vector<char> used(num_verts_, 0);
    
    for (int start = 0; start < num_verts_; ++start) {
      if (used[start]) continue;
      
      // Count face size first
      int size = 1;
      int curr = verts_[start].next;
      while (curr != start) {
        curr = verts_[curr].next;
        size++;
        if (size > num_verts_) break;  // Safety
      }
      
      if (size < 3) {
        used[start] = 1;
        continue;
      }
      
      // Extract face vertices
      std::vector<int> face;
      face.reserve(size);
      curr = start;
      do {
        face.push_back(curr);
        used[curr] = 1;
        curr = verts_[curr].next;
      } while (curr != start);
      
      // Check CCW orientation
      double area = 0;
      for (size_t j = 0; j < face.size(); ++j) {
        size_t k = (j + 1) % face.size();
        area += verts_[face[j]].x * verts_[face[k]].y -
                verts_[face[k]].x * verts_[face[j]].y;
      }
      if (area > kEps) {
        num_faces++;
        triangulate_monotone(face);
      }
    }
  }
  
  void triangulate_monotone(const std::vector<int>& face) {
    const int m = static_cast<int>(face.size());
    
    if (m < 3) return;
    if (m == 3) {
      triangles.push_back({face[0] % n_, face[1] % n_, face[2] % n_});
      return;
    }
    
    // Find top vertex
    int top = 0;
    for (int i = 1; i < m; ++i) {
      if (verts_[face[i]].y > verts_[face[top]].y + kEps ||
          (std::abs(verts_[face[i]].y - verts_[face[top]].y) < kEps &&
           verts_[face[i]].x < verts_[face[top]].x)) {
        top = i;
      }
    }
    
    // Sort by y descending
    std::vector<std::pair<int, int>> sorted(m);
    for (int i = 0; i < m; ++i) {
      sorted[i] = {i, face[i]};
    }
    std::sort(sorted.begin(), sorted.end(), [this](const auto& a, const auto& b) {
      double ya = verts_[a.second].y, yb = verts_[b.second].y;
      if (std::abs(ya - yb) > kEps) return ya > yb;
      return verts_[a.second].x < verts_[b.second].x;
    });
    
    int bot = sorted.back().first;
    
    // Assign sides
    std::vector<int8_t> side(m);
    for (int i = top; ; i = (i + 1) % m) {
      side[i] = 0;
      if (i == bot) break;
    }
    for (int i = top; ; i = (i - 1 + m) % m) {
      side[i] = 1;
      if (i == bot) break;
    }
    
    // Stack-based monotone triangulation
    std::vector<int> stk;
    stk.push_back(sorted[0].first);
    stk.push_back(sorted[1].first);
    
    for (int i = 2; i < m; ++i) {
      int vi = sorted[i].first;
      
      if (side[vi] != side[stk.back()]) {
        while (stk.size() > 1) {
          int a = stk[stk.size() - 2];
          int b = stk[stk.size() - 1];
          triangles.push_back({face[a] % n_, face[b] % n_, face[vi] % n_});
          stk.pop_back();
        }
        stk.pop_back();
        stk.push_back(sorted[i - 1].first);
        stk.push_back(vi);
      } else {
        int last = stk.back();
        stk.pop_back();
        
        while (!stk.empty()) {
          int a = stk.back();
          double c = (verts_[face[a]].x - verts_[face[vi]].x) *
                     (verts_[face[last]].y - verts_[face[vi]].y) -
                     (verts_[face[a]].y - verts_[face[vi]].y) *
                     (verts_[face[last]].x - verts_[face[vi]].x);
          bool ok = (side[vi] == 0) ? (c > kEps) : (c < -kEps);
          if (!ok) break;
          triangles.push_back({face[a] % n_, face[last] % n_, face[vi] % n_});
          last = stk.back();
          stk.pop_back();
        }
        stk.push_back(last);
        stk.push_back(vi);
      }
    }
  }
};

}  // namespace fast_linked
