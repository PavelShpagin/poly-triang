#pragma once
/**
 * Reflex Triangulation Algorithm with Spatial Hashing
 * 
 * Complexity: O(n + r log r) expected for well-distributed polygons
 * 
 * Key optimizations:
 * 1. Convex fast-path: O(n) fan triangulation when r=0
 * 2. Spatial hashing: Only check nearby reflex vertices for point-in-triangle
 * 3. Reflex set maintenance: Track which vertices are still reflex
 */

#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

namespace reflex {

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

inline double cross(const Point& p0, const Point& p1, const Point& p2) {
    return (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x);
}

inline double signed_area2(const Point& a, const Point& b, const Point& c) {
    return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
}

inline bool point_in_triangle_strict(const Point& p, const Point& a, const Point& b, const Point& c) {
    double d1 = signed_area2(p, a, b);
    double d2 = signed_area2(p, b, c);
    double d3 = signed_area2(p, c, a);
    
    bool has_neg = (d1 < -1e-10) || (d2 < -1e-10) || (d3 < -1e-10);
    bool has_pos = (d1 > 1e-10) || (d2 > 1e-10) || (d3 > 1e-10);
    
    return !(has_neg && has_pos);
}

// Simple spatial hash for fast point queries
class SpatialHash {
    double cell_size;
    double inv_cell_size;
    std::unordered_map<int64_t, std::vector<int>> cells;
    
    int64_t hash_key(double x, double y) const {
        int64_t ix = static_cast<int64_t>(std::floor(x * inv_cell_size));
        int64_t iy = static_cast<int64_t>(std::floor(y * inv_cell_size));
        return (ix * 73856093) ^ (iy * 19349663);  // Hash combine
    }
    
public:
    SpatialHash(double cell_sz = 1.0) : cell_size(cell_sz), inv_cell_size(1.0 / cell_sz) {}
    
    void clear() { cells.clear(); }
    
    void insert(int idx, double x, double y) {
        cells[hash_key(x, y)].push_back(idx);
    }
    
    void remove(int idx, double x, double y) {
        auto& cell = cells[hash_key(x, y)];
        cell.erase(std::remove(cell.begin(), cell.end(), idx), cell.end());
    }
    
    // Query all points in cells overlapping the bounding box
    template<typename Func>
    void query_bbox(double minx, double miny, double maxx, double maxy, Func&& func) const {
        int64_t ix0 = static_cast<int64_t>(std::floor(minx * inv_cell_size));
        int64_t iy0 = static_cast<int64_t>(std::floor(miny * inv_cell_size));
        int64_t ix1 = static_cast<int64_t>(std::floor(maxx * inv_cell_size));
        int64_t iy1 = static_cast<int64_t>(std::floor(maxy * inv_cell_size));
        
        for (int64_t ix = ix0; ix <= ix1; ix++) {
            for (int64_t iy = iy0; iy <= iy1; iy++) {
                int64_t key = (ix * 73856093) ^ (iy * 19349663);
                auto it = cells.find(key);
                if (it != cells.end()) {
                    for (int idx : it->second) {
                        func(idx);
                    }
                }
            }
        }
    }
};

class ReflexTriangulator {
public:
    std::vector<Triangle> triangles;
    int num_reflex = 0;
    
    void triangulate(const std::vector<Point>& input_polygon) {
        triangles.clear();
        num_reflex = 0;
        
        int n = input_polygon.size();
        if (n < 3) return;
        
        // Copy polygon
        std::vector<Point> poly = input_polygon;
        for (int i = 0; i < n; i++) {
            poly[i].index = i;
        }
        
        // Build doubly linked list
        std::vector<int> prev(n), next(n);
        for (int i = 0; i < n; i++) {
            prev[i] = (i - 1 + n) % n;
            next[i] = (i + 1) % n;
        }
        
        // Compute bounding box and cell size
        double minx = poly[0].x, maxx = poly[0].x;
        double miny = poly[0].y, maxy = poly[0].y;
        for (int i = 1; i < n; i++) {
            minx = std::min(minx, poly[i].x);
            maxx = std::max(maxx, poly[i].x);
            miny = std::min(miny, poly[i].y);
            maxy = std::max(maxy, poly[i].y);
        }
        double diag = std::sqrt((maxx - minx) * (maxx - minx) + (maxy - miny) * (maxy - miny));
        double cell_size = std::max(diag / std::sqrt(n), 1e-6);
        
        // Find reflex vertices and build spatial hash
        std::unordered_set<int> reflex_set;
        SpatialHash spatial_hash(cell_size);
        
        for (int i = 0; i < n; i++) {
            if (cross(poly[prev[i]], poly[i], poly[next[i]]) < -1e-10) {
                reflex_set.insert(i);
                spatial_hash.insert(i, poly[i].x, poly[i].y);
            }
        }
        num_reflex = reflex_set.size();
        
        // Special case: convex polygon - fan triangulation O(n)
        if (reflex_set.empty()) {
            for (int i = 1; i < n - 1; i++) {
                triangles.emplace_back(0, i, i + 1);
            }
            return;
        }
        
        // Ear clipping with spatial hash optimization
        auto is_convex = [&](int i) -> bool {
            return cross(poly[prev[i]], poly[i], poly[next[i]]) > 1e-10;
        };
        
        auto is_ear = [&](int i) -> bool {
            if (!is_convex(i)) return false;
            
            const Point& a = poly[prev[i]];
            const Point& b = poly[i];
            const Point& c = poly[next[i]];
            
            // Compute bounding box of triangle
            double tri_minx = std::min({a.x, b.x, c.x});
            double tri_maxx = std::max({a.x, b.x, c.x});
            double tri_miny = std::min({a.y, b.y, c.y});
            double tri_maxy = std::max({a.y, b.y, c.y});
            
            // Query spatial hash for nearby reflex vertices
            bool has_interior_point = false;
            spatial_hash.query_bbox(tri_minx, tri_miny, tri_maxx, tri_maxy, [&](int r) {
                if (has_interior_point) return;
                if (r == prev[i] || r == i || r == next[i]) return;
                if (!reflex_set.count(r)) return;  // Already removed
                if (point_in_triangle_strict(poly[r], a, b, c)) {
                    has_interior_point = true;
                }
            });
            
            return !has_interior_point;
        };
        
        int remaining = n;
        int start = 0;
        int max_iter = n * n;
        int iter = 0;
        
        while (remaining > 3 && iter < max_iter) {
            iter++;
            int current = start;
            bool found = false;
            
            do {
                if (is_ear(current)) {
                    triangles.emplace_back(prev[current], current, next[current]);
                    
                    int p = prev[current];
                    int nx = next[current];
                    
                    next[p] = nx;
                    prev[nx] = p;
                    
                    // Remove from reflex set and spatial hash
                    if (reflex_set.count(current)) {
                        reflex_set.erase(current);
                        spatial_hash.remove(current, poly[current].x, poly[current].y);
                    }
                    
                    // Update reflex status of neighbors
                    if (reflex_set.count(p) && cross(poly[prev[p]], poly[p], poly[next[p]]) > 1e-10) {
                        reflex_set.erase(p);
                        spatial_hash.remove(p, poly[p].x, poly[p].y);
                    }
                    if (reflex_set.count(nx) && cross(poly[prev[nx]], poly[nx], poly[next[nx]]) > 1e-10) {
                        reflex_set.erase(nx);
                        spatial_hash.remove(nx, poly[nx].x, poly[nx].y);
                    }
                    
                    remaining--;
                    start = nx;
                    found = true;
                    break;
                }
                current = next[current];
            } while (current != start);
            
            if (!found) {
                // Fallback: find any convex vertex
                current = start;
                do {
                    if (is_convex(current)) {
                        triangles.emplace_back(prev[current], current, next[current]);
                        int p = prev[current];
                        int nx = next[current];
                        next[p] = nx;
                        prev[nx] = p;
                        if (reflex_set.count(current)) {
                            reflex_set.erase(current);
                            spatial_hash.remove(current, poly[current].x, poly[current].y);
                        }
                        remaining--;
                        start = nx;
                        found = true;
                        break;
                    }
                    current = next[current];
                } while (current != start);
                
                if (!found) {
                    // Force progress
                    triangles.emplace_back(prev[start], start, next[start]);
                    int p = prev[start];
                    int nx = next[start];
                    next[p] = nx;
                    prev[nx] = p;
                    if (reflex_set.count(start)) {
                        reflex_set.erase(start);
                        spatial_hash.remove(start, poly[start].x, poly[start].y);
                    }
                    remaining--;
                    start = nx;
                }
            }
        }
        
        if (remaining == 3) {
            int v = start;
            triangles.emplace_back(prev[v], v, next[v]);
        }
    }
};

} // namespace reflex
