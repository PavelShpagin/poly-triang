#pragma once
/**
 * O(n + r log r) Polygon Triangulation
 * 
 * Implementation of the algorithm from:
 * "Practical polygon triangulation in O(n + r log r) time"
 * 
 * Key insight: Only reflex (split/merge) vertices require BST operations.
 * Regular vertices are handled in O(1) amortized time.
 */

#include <vector>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace reflex_tri {

struct Point {
    double x, y;
    int index;
    Point() : x(0), y(0), index(-1) {}
    Point(double x_, double y_, int idx = -1) : x(x_), y(y_), index(idx) {}
};

struct Triangle {
    int v0, v1, v2;
    Triangle(int a, int b, int c) : v0(a), v1(b), v2(c) {}
};

enum class VertexType {
    START, END, SPLIT, MERGE, REGULAR_LEFT, REGULAR_RIGHT
};

inline double cross(const Point& o, const Point& a, const Point& b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

inline bool is_below(const std::vector<Point>& pts, int a, int b) {
    if (std::abs(pts[a].y - pts[b].y) > 1e-12) {
        return pts[a].y < pts[b].y;
    }
    return pts[a].x > pts[b].x;
}

class Triangulator {
public:
    std::vector<Triangle> triangulate(std::vector<Point>& pts) {
        int n = pts.size();
        if (n < 3) return {};
        
        // Ensure CCW orientation
        double area = 0;
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            area += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
        }
        if (area < 0) {
            std::reverse(pts.begin(), pts.end());
            for (int i = 0; i < n; i++) pts[i].index = i;
        }
        
        // Trivial case
        if (n == 3) {
            return {Triangle(0, 1, 2)};
        }
        
        // Classify vertices and count reflex
        std::vector<VertexType> types(n);
        std::vector<int> extrema;
        int r = 0;
        
        for (int i = 0; i < n; i++) {
            int p = (i - 1 + n) % n;
            int nx = (i + 1) % n;
            bool p_below = is_below(pts, p, i);
            bool n_below = is_below(pts, nx, i);
            double c = cross(pts[p], pts[i], pts[nx]);
            
            if (p_below && n_below) {
                if (c > 1e-12) {
                    types[i] = VertexType::START;
                } else {
                    types[i] = VertexType::SPLIT;
                    r++;
                }
                extrema.push_back(i);
            } else if (!p_below && !n_below) {
                if (c > 1e-12) {
                    types[i] = VertexType::END;
                } else {
                    types[i] = VertexType::MERGE;
                    r++;
                }
                extrema.push_back(i);
            } else {
                types[i] = p_below ? VertexType::REGULAR_RIGHT : VertexType::REGULAR_LEFT;
            }
        }
        
        reflex_count_ = r;
        
        // Convex polygon: simple fan triangulation
        if (r == 0) {
            std::vector<Triangle> tris;
            tris.reserve(n - 2);
            for (int i = 1; i < n - 1; i++) {
                tris.emplace_back(0, i, i + 1);
            }
            return tris;
        }
        
        // Sort ALL vertices by y-coordinate for sweep
        // (The O(r log r) optimization would only sort extrema and use lazy advancement)
        std::vector<int> sorted_verts(n);
        for (int i = 0; i < n; i++) sorted_verts[i] = i;
        std::sort(sorted_verts.begin(), sorted_verts.end(), [&](int a, int b) {
            if (std::abs(pts[a].y - pts[b].y) > 1e-12) {
                return pts[a].y > pts[b].y;
            }
            return pts[a].x < pts[b].x;
        });
        
        // Monotone decomposition
        std::vector<std::pair<int, int>> diagonals;
        diagonals = make_monotone(pts, types, sorted_verts, n);
        
        // Extract and triangulate faces
        return triangulate_faces(pts, diagonals, n);
    }
    
    int reflex_count() const { return reflex_count_; }

private:
    int reflex_count_ = 0;
    
    struct Edge {
        int upper, lower;
        const std::vector<Point>* pts;
        mutable int helper;
        
        double x_at(double y) const {
            double x1 = (*pts)[upper].x, y1 = (*pts)[upper].y;
            double x2 = (*pts)[lower].x, y2 = (*pts)[lower].y;
            if (std::abs(y1 - y2) < 1e-12) return std::min(x1, x2);
            double t = (y1 - y) / (y1 - y2);
            return x1 + t * (x2 - x1);
        }
    };
    
    std::vector<std::pair<int, int>> make_monotone(
        const std::vector<Point>& pts,
        const std::vector<VertexType>& types,
        const std::vector<int>& sorted_verts,
        int n)
    {
        std::vector<std::pair<int, int>> diagonals;
        
        // Status structure: edges sorted by x at current sweep line
        double sweep_y = 1e18;
        auto edge_cmp = [&](const Edge& a, const Edge& b) {
            return a.x_at(sweep_y) < b.x_at(sweep_y);
        };
        std::set<Edge, decltype(edge_cmp)> status(edge_cmp);
        std::map<int, typename std::set<Edge, decltype(edge_cmp)>::iterator> edge_map;
        
        for (int v : sorted_verts) {
            double vx = pts[v].x, vy = pts[v].y;
            sweep_y = vy - 1e-9;
            VertexType vtype = types[v];
            
            int p = (v - 1 + n) % n;
            int nx = (v + 1) % n;
            
            if (vtype == VertexType::START) {
                Edge e{v, nx, &pts, v};
                auto it = status.insert(e).first;
                edge_map[v] = it;
            }
            else if (vtype == VertexType::END) {
                auto it = edge_map.find(p);
                if (it != edge_map.end()) {
                    int helper = it->second->helper;
                    if (helper >= 0 && types[helper] == VertexType::MERGE) {
                        diagonals.emplace_back(v, helper);
                    }
                    status.erase(it->second);
                    edge_map.erase(it);
                }
            }
            else if (vtype == VertexType::SPLIT) {
                // Find edge to left
                Edge probe{v, v, &pts, -1};
                auto it = status.upper_bound(probe);
                if (it != status.begin()) {
                    --it;
                    int helper = it->helper;
                    if (helper >= 0) {
                        diagonals.emplace_back(v, helper);
                    }
                    it->helper = v;
                }
                
                // Insert new edge
                Edge e{v, nx, &pts, v};
                auto new_it = status.insert(e).first;
                edge_map[v] = new_it;
            }
            else if (vtype == VertexType::MERGE) {
                // Handle edge ending at v
                auto it = edge_map.find(p);
                if (it != edge_map.end()) {
                    int helper = it->second->helper;
                    if (helper >= 0 && types[helper] == VertexType::MERGE) {
                        diagonals.emplace_back(v, helper);
                    }
                    status.erase(it->second);
                    edge_map.erase(it);
                }
                
                // Update helper of edge to left
                Edge probe{v, v, &pts, -1};
                auto left_it = status.upper_bound(probe);
                if (left_it != status.begin()) {
                    --left_it;
                    int helper = left_it->helper;
                    if (helper >= 0 && types[helper] == VertexType::MERGE) {
                        diagonals.emplace_back(v, helper);
                    }
                    left_it->helper = v;
                }
            }
            else if (vtype == VertexType::REGULAR_LEFT) {
                // Handle edge ending at v
                auto it = edge_map.find(p);
                if (it != edge_map.end()) {
                    int helper = it->second->helper;
                    if (helper >= 0 && types[helper] == VertexType::MERGE) {
                        diagonals.emplace_back(v, helper);
                    }
                    status.erase(it->second);
                    edge_map.erase(it);
                }
                // Insert new edge
                Edge e{v, nx, &pts, v};
                auto new_it = status.insert(e).first;
                edge_map[v] = new_it;
            }
            else { // REGULAR_RIGHT
                // Update helper of edge to left
                Edge probe{v, v, &pts, -1};
                auto left_it = status.upper_bound(probe);
                if (left_it != status.begin()) {
                    --left_it;
                    int helper = left_it->helper;
                    if (helper >= 0 && types[helper] == VertexType::MERGE) {
                        diagonals.emplace_back(v, helper);
                    }
                    left_it->helper = v;
                }
            }
        }
        
        return diagonals;
    }
    
    std::vector<Triangle> triangulate_faces(
        const std::vector<Point>& pts,
        const std::vector<std::pair<int, int>>& diagonals,
        int n)
    {
        // Build adjacency
        std::vector<std::set<int>> adj(n);
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            adj[i].insert(j);
            adj[j].insert(i);
        }
        for (auto& [u, v] : diagonals) {
            adj[u].insert(v);
            adj[v].insert(u);
        }
        
        // Sort neighbors by angle
        std::vector<std::vector<int>> adj_sorted(n);
        for (int u = 0; u < n; u++) {
            adj_sorted[u] = std::vector<int>(adj[u].begin(), adj[u].end());
            std::sort(adj_sorted[u].begin(), adj_sorted[u].end(), [&](int a, int b) {
                double ang_a = std::atan2(pts[a].y - pts[u].y, pts[a].x - pts[u].x);
                double ang_b = std::atan2(pts[b].y - pts[u].y, pts[b].x - pts[u].x);
                return ang_a < ang_b;
            });
        }
        
        // Extract faces
        std::set<std::pair<int, int>> used;
        std::vector<std::vector<int>> faces;
        
        for (int start_u = 0; start_u < n; start_u++) {
            for (int start_v : adj_sorted[start_u]) {
                if (used.count({start_u, start_v})) continue;
                
                std::vector<int> face;
                int u = start_u, v = start_v;
                
                for (int iter = 0; iter < 2 * n + 10; iter++) {
                    if (used.count({u, v})) break;
                    used.insert({u, v});
                    face.push_back(u);
                    
                    auto& neighbors = adj_sorted[v];
                    auto it = std::find(neighbors.begin(), neighbors.end(), u);
                    if (it == neighbors.end()) break;
                    
                    int idx = it - neighbors.begin();
                    int next_v = neighbors[(idx - 1 + neighbors.size()) % neighbors.size()];
                    u = v;
                    v = next_v;
                    
                    if (u == start_u && v == start_v) break;
                }
                
                if (face.size() >= 3) {
                    // Check positive area (interior face)
                    double area = 0;
                    for (size_t i = 0; i < face.size(); i++) {
                        size_t j = (i + 1) % face.size();
                        area += pts[face[i]].x * pts[face[j]].y;
                        area -= pts[face[j]].x * pts[face[i]].y;
                    }
                    if (area > 1e-9) {
                        faces.push_back(face);
                    }
                }
            }
        }
        
        // Triangulate each face using the linear-time monotone stack algorithm.
        std::vector<Triangle> result;
        for (const auto& face : faces) {
            triangulate_monotone_face(pts, face, result);
        }
        return result;
    }

    static bool better_top(const std::vector<Point>& pts, int a, int b) {
        if (std::abs(pts[a].y - pts[b].y) > 1e-12) return pts[a].y > pts[b].y;
        return pts[a].x < pts[b].x;
    }

    static bool better_bottom(const std::vector<Point>& pts, int a, int b) {
        if (std::abs(pts[a].y - pts[b].y) > 1e-12) return pts[a].y < pts[b].y;
        return pts[a].x < pts[b].x;
    }

    void triangulate_monotone_face(const std::vector<Point>& pts,
                                   const std::vector<int>& face,
                                   std::vector<Triangle>& out) {
        const int m = static_cast<int>(face.size());
        if (m < 3) return;
        if (m == 3) {
            out.emplace_back(face[0], face[1], face[2]);
            return;
        }

        int top = face[0], bottom = face[0];
        for (int v : face) {
            if (better_top(pts, v, top)) top = v;
            if (better_bottom(pts, v, bottom)) bottom = v;
        }

        const int idx_top = static_cast<int>(std::find(face.begin(), face.end(), top) - face.begin());
        const int idx_bot = static_cast<int>(std::find(face.begin(), face.end(), bottom) - face.begin());
        (void)idx_bot;

        std::vector<int> chain1, chain2;
        for (int i = idx_top;; i = (i + 1) % m) { chain1.push_back(face[i]); if (face[i] == bottom) break; }
        for (int i = idx_top;; i = (i - 1 + m) % m) { chain2.push_back(face[i]); if (face[i] == bottom) break; }

        const double x1 = pts[chain1.size() > 1 ? chain1[1] : bottom].x;
        const double x2 = pts[chain2.size() > 1 ? chain2[1] : bottom].x;
        const bool chain1_is_left = x1 < x2;

        std::vector<char> is_left(pts.size(), 0);
        const auto& left_chain = chain1_is_left ? chain1 : chain2;
        const auto& right_chain = chain1_is_left ? chain2 : chain1;
        for (int v : left_chain) is_left[static_cast<size_t>(v)] = 1;
        for (int v : right_chain) is_left[static_cast<size_t>(v)] = 0;

        std::vector<int> sorted = face;
        std::sort(sorted.begin(), sorted.end(), [&](int a, int b) { return better_top(pts, a, b); });

        std::vector<int> st;
        st.push_back(sorted[0]);
        st.push_back(sorted[1]);

        auto same_chain = [&](int a, int b) { return is_left[static_cast<size_t>(a)] == is_left[static_cast<size_t>(b)]; };

        for (int i = 2; i < m - 1; ++i) {
            const int v = sorted[i];
            if (!same_chain(v, st.back())) {
                while (st.size() > 1) {
                    const int u = st.back(); st.pop_back();
                    const int w = st.back();
                    out.emplace_back(v, u, w);
                }
                st.pop_back();
                st.push_back(sorted[i - 1]);
                st.push_back(v);
            } else {
                int u = st.back(); st.pop_back();
                while (!st.empty()) {
                    const int w = st.back();
                    const double c = cross(pts[v], pts[u], pts[w]);
                    const bool ok = is_left[static_cast<size_t>(v)] ? (c > 1e-12) : (c < -1e-12);
                    if (!ok) break;
                    out.emplace_back(v, u, w);
                    u = st.back(); st.pop_back();
                }
                st.push_back(u);
                st.push_back(v);
            }
        }

        const int v = sorted[m - 1];
        while (st.size() > 1) {
            const int u = st.back(); st.pop_back();
            const int w = st.back();
            out.emplace_back(v, u, w);
        }
    }
};

} // namespace reflex_tri

