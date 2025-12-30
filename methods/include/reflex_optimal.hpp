#pragma once
/**
 * Optimized Polygon Triangulation via Monotone Decomposition
 * 
 * Complexity:
 * - O(n) for convex polygons (fast path)
 * - O(n log n) for general polygons (optimized edge-based sweep)
 * 
 * Optimizations:
 * - Binary search for edge lookup O(log n) per query
 * - Cross-product angle comparison (no trig functions)
 * - Cache-friendly data structures
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <unordered_map>

namespace reflex_opt {

struct Point {
  double x = 0.0;
  double y = 0.0;
  int index = -1;
};

struct Triangle { int v0, v1, v2; };

namespace detail {
constexpr double kEps = 1e-9;

inline double cross(const Point& o, const Point& a, const Point& b) {
  return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

inline double signed_area(const std::vector<Point>& pts) {
  double a = 0.0;
  const int n = static_cast<int>(pts.size());
  for (int i = 0; i < n; ++i) {
    int j = (i + 1) % n;
    a += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
  }
  return 0.5 * a;
}

inline bool below(const Point& a, const Point& b) {
  return (a.y < b.y - kEps) || (std::abs(a.y - b.y) < kEps && a.x > b.x);
}

}  // namespace detail

class Triangulator {
public:
  std::vector<Triangle> triangulate(std::vector<Point>& pts) {
    triangles_.clear();
    diagonals_.clear();
    reflex_count_ = 0;
    face_count_ = 0;

    n_ = static_cast<int>(pts.size());
    if (n_ < 3) return {};

    for (int i = 0; i < n_; ++i) pts[i].index = i;
    
    if (detail::signed_area(pts) < 0.0) {
      std::reverse(pts.begin(), pts.end());
    }

    pts_ = &pts;

    for (int i = 0; i < n_; ++i) {
      if (is_reflex(i)) ++reflex_count_;
    }

    if (n_ == 3) {
      triangles_.push_back({pts[0].index, pts[1].index, pts[2].index});
      face_count_ = 1;
      return triangles_;
    }

    // Fast path for convex - O(n)
    if (reflex_count_ == 0) {
      triangles_.reserve(n_ - 2);
      for (int i = 1; i < n_ - 1; ++i) {
        triangles_.push_back({pts[0].index, pts[i].index, pts[i + 1].index});
      }
      face_count_ = 1;
      return triangles_;
    }

    classify_vertices();
    decompose();
    triangulate_faces();

    for (auto& tri : triangles_) {
      tri.v0 = pts[tri.v0].index;
      tri.v1 = pts[tri.v1].index;
      tri.v2 = pts[tri.v2].index;
    }

    return triangles_;
  }

  int reflex_count() const { return reflex_count_; }
  int diagonal_count() const { return static_cast<int>(diagonals_.size()); }
  int face_count() const { return face_count_; }

private:
  enum VType : uint8_t { START, END, SPLIT, MERGE, REGULAR };
  
  const std::vector<Point>* pts_ = nullptr;
  int n_ = 0;
  int reflex_count_ = 0;
  int face_count_ = 0;
  std::vector<VType> types_;
  std::vector<std::pair<int, int>> diagonals_;
  std::vector<Triangle> triangles_;

  int prev(int i) const { return (i + n_ - 1) % n_; }
  int next(int i) const { return (i + 1) % n_; }
  const Point& pt(int i) const { return (*pts_)[i]; }
  
  bool is_reflex(int i) const {
    return detail::cross(pt(prev(i)), pt(i), pt(next(i))) < -detail::kEps;
  }

  void classify_vertices() {
    types_.resize(n_);
    for (int i = 0; i < n_; ++i) {
      bool prev_below = detail::below(pt(prev(i)), pt(i));
      bool next_below = detail::below(pt(next(i)), pt(i));
      
      if (prev_below && next_below) {
        types_[i] = is_reflex(i) ? SPLIT : START;
      } else if (!prev_below && !next_below) {
        types_[i] = is_reflex(i) ? MERGE : END;
      } else {
        types_[i] = REGULAR;
      }
    }
  }

  std::unordered_map<int, int> edge_helper_;
  
  double edge_x_at(int v, double y) const {
    const Point& p = pt(v);
    const Point& q = pt(next(v));
    double dy = p.y - q.y;
    if (std::abs(dy) < detail::kEps) return std::min(p.x, q.x);
    return p.x + (p.y - y) / dy * (q.x - p.x);
  }

  void decompose() {
    edge_helper_.clear();
    edge_helper_.reserve(n_);
    
    std::vector<int> sorted(n_);
    for (int i = 0; i < n_; ++i) sorted[i] = i;
    std::sort(sorted.begin(), sorted.end(), [this](int a, int b) {
      double ya = pt(a).y, yb = pt(b).y;
      if (std::abs(ya - yb) > detail::kEps) return ya > yb;
      return pt(a).x < pt(b).x;
    });
    
    std::vector<int> active;
    active.reserve(n_ / 4);
    double sweep_y;
    
    auto get_x = [&](int e) { return edge_x_at(e, sweep_y); };
    
    auto find_left = [&](double x) -> int {
      if (active.empty()) return -1;
      int lo = 0, hi = static_cast<int>(active.size()) - 1, result = -1;
      while (lo <= hi) {
        int mid = (lo + hi) / 2;
        if (get_x(active[mid]) < x - detail::kEps) {
          result = active[mid];
          lo = mid + 1;
        } else {
          hi = mid - 1;
        }
      }
      return result;
    };
    
    auto insert_edge = [&](int e) {
      double ex = get_x(e);
      auto it = std::lower_bound(active.begin(), active.end(), ex,
        [&](int edge, double x) { return get_x(edge) < x; });
      active.insert(it, e);
    };
    
    auto remove_edge = [&](int e) {
      double ex = get_x(e);
      auto it = std::lower_bound(active.begin(), active.end(), ex,
        [&](int edge, double x) { return get_x(edge) < x; });
      while (it != active.end() && *it != e) ++it;
      if (it != active.end()) active.erase(it);
    };
    
    for (int v : sorted) {
      sweep_y = pt(v).y;
      
      switch (types_[v]) {
        case START: {
          insert_edge(v);
          edge_helper_[v] = v;
          break;
        }
        case END: {
          int e = prev(v);
          auto it = edge_helper_.find(e);
          if (it != edge_helper_.end() && types_[it->second] == MERGE) {
            add_diagonal(v, it->second);
          }
          remove_edge(e);
          edge_helper_.erase(e);
          break;
        }
        case SPLIT: {
          int left = find_left(pt(v).x);
          if (left >= 0) {
            auto it = edge_helper_.find(left);
            if (it != edge_helper_.end()) {
              add_diagonal(v, it->second);
              it->second = v;
            }
          }
          insert_edge(v);
          edge_helper_[v] = v;
          break;
        }
        case MERGE: {
          int e = prev(v);
          auto it = edge_helper_.find(e);
          if (it != edge_helper_.end() && types_[it->second] == MERGE) {
            add_diagonal(v, it->second);
          }
          remove_edge(e);
          edge_helper_.erase(e);
          int left = find_left(pt(v).x);
          if (left >= 0) {
            auto it2 = edge_helper_.find(left);
            if (it2 != edge_helper_.end()) {
              if (types_[it2->second] == MERGE) {
                add_diagonal(v, it2->second);
              }
              it2->second = v;
            }
          }
          break;
        }
        case REGULAR: {
          bool interior_right = pt(prev(v)).y > pt(v).y + detail::kEps ||
            (std::abs(pt(prev(v)).y - pt(v).y) < detail::kEps && pt(prev(v)).x < pt(v).x);
          
          if (interior_right) {
            int e = prev(v);
            auto it = edge_helper_.find(e);
            if (it != edge_helper_.end() && types_[it->second] == MERGE) {
              add_diagonal(v, it->second);
            }
            remove_edge(e);
            edge_helper_.erase(e);
            insert_edge(v);
            edge_helper_[v] = v;
          } else {
            int left = find_left(pt(v).x);
            if (left >= 0) {
              auto it = edge_helper_.find(left);
              if (it != edge_helper_.end()) {
                if (types_[it->second] == MERGE) {
                  add_diagonal(v, it->second);
                }
                it->second = v;
              }
            }
          }
          break;
        }
      }
    }
    
    for (auto& d : diagonals_) {
      if (d.first > d.second) std::swap(d.first, d.second);
    }
    std::sort(diagonals_.begin(), diagonals_.end());
    diagonals_.erase(std::unique(diagonals_.begin(), diagonals_.end()), diagonals_.end());
  }
  
  void add_diagonal(int u, int v) {
    if (u == v) return;
    int diff = std::abs(u - v);
    if (diff == 1 || diff == n_ - 1) return;
    diagonals_.push_back({u, v});
  }

  static int angle_cmp(double ux, double uy, double ax, double ay, double bx, double by) {
    double dax = ax - ux, day = ay - uy;
    double dbx = bx - ux, dby = by - uy;
    int qa = (day < 0 || (day == 0 && dax < 0)) ? 1 : 0;
    int qb = (dby < 0 || (dby == 0 && dbx < 0)) ? 1 : 0;
    if (qa != qb) return qa - qb;
    double cross = dax * dby - day * dbx;
    if (cross > detail::kEps) return -1;
    if (cross < -detail::kEps) return 1;
    return 0;
  }

  void triangulate_faces() {
    const auto& pts = *pts_;
    
    std::vector<std::vector<int>> adj(n_);
    for (int i = 0; i < n_; ++i) {
      adj[i].reserve(4);
      adj[i].push_back(next(i));
      adj[i].push_back(prev(i));
    }
    for (const auto& d : diagonals_) {
      adj[d.first].push_back(d.second);
      adj[d.second].push_back(d.first);
    }
    
    for (int u = 0; u < n_; ++u) {
      double ux = pts[u].x, uy = pts[u].y;
      std::sort(adj[u].begin(), adj[u].end(), [&](int a, int b) {
        return angle_cmp(ux, uy, pts[a].x, pts[a].y, pts[b].x, pts[b].y) < 0;
      });
    }
    
    std::vector<std::vector<bool>> used(n_);
    for (int i = 0; i < n_; ++i) {
      used[i].assign(adj[i].size(), false);
    }
    
    auto find_idx = [&](int u, int v) -> int {
      const auto& nb = adj[u];
      for (int i = 0; i < static_cast<int>(nb.size()); ++i) {
        if (nb[i] == v) return i;
      }
      return -1;
    };
    
    std::vector<std::vector<int>> faces;
    faces.reserve(diagonals_.size() + 1);
    
    for (int u = 0; u < n_; ++u) {
      for (int vi = 0; vi < static_cast<int>(adj[u].size()); ++vi) {
        if (used[u][vi]) continue;
        
        std::vector<int> face;
        face.reserve(8);
        int cu = u, cvi = vi;
        
        while (!used[cu][cvi]) {
          used[cu][cvi] = true;
          face.push_back(cu);
          
          int cv = adj[cu][cvi];
          int back_idx = find_idx(cv, cu);
          if (back_idx < 0) break;
          
          int next_idx = (back_idx + static_cast<int>(adj[cv].size()) - 1) % static_cast<int>(adj[cv].size());
          cu = cv;
          cvi = next_idx;
          
          if (cu == u && adj[cu][cvi] == adj[u][vi]) break;
        }
        
        if (face.size() >= 3) {
          double area = 0;
          for (size_t i = 0; i < face.size(); ++i) {
            size_t j = (i + 1) % face.size();
            area += pts[face[i]].x * pts[face[j]].y - pts[face[j]].x * pts[face[i]].y;
          }
          if (area > detail::kEps) {
            faces.push_back(std::move(face));
          }
        }
      }
    }
    
    face_count_ = static_cast<int>(faces.size());
    triangles_.reserve(n_ - 2);
    
    for (const auto& face : faces) {
      triangulate_monotone(face);
    }
  }

  void triangulate_monotone(const std::vector<int>& face) {
    const auto& pts = *pts_;
    const int m = static_cast<int>(face.size());
    
    if (m < 3) return;
    if (m == 3) {
      triangles_.push_back({face[0], face[1], face[2]});
      return;
    }
    
    int top_idx = 0;
    for (int i = 1; i < m; ++i) {
      if (pts[face[i]].y > pts[face[top_idx]].y + detail::kEps ||
          (std::abs(pts[face[i]].y - pts[face[top_idx]].y) < detail::kEps &&
           pts[face[i]].x < pts[face[top_idx]].x)) {
        top_idx = i;
      }
    }
    
    std::vector<std::pair<int, int>> sorted(m);
    for (int i = 0; i < m; ++i) sorted[i] = {i, face[i]};
    std::sort(sorted.begin(), sorted.end(), [&](const auto& a, const auto& b) {
      double ya = pts[a.second].y, yb = pts[b.second].y;
      if (std::abs(ya - yb) > detail::kEps) return ya > yb;
      return pts[a.second].x < pts[b.second].x;
    });
    
    int bot_idx = sorted.back().first;
    
    std::vector<int8_t> side(m);
    for (int i = top_idx; ; i = (i + 1) % m) {
      side[i] = 0;
      if (i == bot_idx) break;
    }
    for (int i = top_idx; ; i = (i + m - 1) % m) {
      side[i] = 1;
      if (i == bot_idx) break;
    }
    side[top_idx] = side[bot_idx] = 0;
    
    std::vector<int> stk;
    stk.reserve(m);
    stk.push_back(sorted[0].first);
    stk.push_back(sorted[1].first);
    
    for (int i = 2; i < m; ++i) {
      int vi = sorted[i].first;
      int v = face[vi];
      
      if (side[vi] != side[stk.back()] || vi == bot_idx) {
        while (stk.size() > 1) {
          int ui = stk.back(); stk.pop_back();
          int wi = stk.back();
          double c = detail::cross(pts[v], pts[face[ui]], pts[face[wi]]);
          if (c > detail::kEps) triangles_.push_back({v, face[ui], face[wi]});
          else if (c < -detail::kEps) triangles_.push_back({v, face[wi], face[ui]});
        }
        stk.pop_back();
        stk.push_back(sorted[i-1].first);
        stk.push_back(vi);
      } else {
        int ui = stk.back(); stk.pop_back();
        while (!stk.empty()) {
          int wi = stk.back();
          double c = detail::cross(pts[v], pts[face[ui]], pts[face[wi]]);
          bool inside = (side[vi] == 0) ? (c < -detail::kEps) : (c > detail::kEps);
          if (!inside) break;
          if (side[vi] == 0) triangles_.push_back({v, face[wi], face[ui]});
          else triangles_.push_back({v, face[ui], face[wi]});
          ui = wi;
          stk.pop_back();
        }
        stk.push_back(ui);
        stk.push_back(vi);
      }
    }
  }
};

}  // namespace reflex_opt
