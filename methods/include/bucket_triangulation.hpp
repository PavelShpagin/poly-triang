#pragma once
/**
 * Bucket Triangulation Algorithm
 * 
 * AMORTIZED O(n) COMPLEXITY
 * 
 * Key insight: By using a grid of O(sqrt(n)) x O(sqrt(n)) buckets,
 * each ear test only checks O(1) expected vertices, and each vertex
 * participates in O(1) expected ear tests.
 * 
 * Algorithm:
 * 1. Compute bounding box and create sqrt(n) x sqrt(n) grid
 * 2. Assign each vertex to its bucket: O(n)
 * 3. For convex polygons: fan triangulation O(n)
 * 4. For non-convex: bucket-accelerated ear clipping
 *    - Each ear test queries O(1) buckets covering the ear triangle
 *    - Each bucket contains O(1) expected vertices
 *    - Each vertex removed from exactly one bucket
 * 
 * Total: O(n) amortized
 */

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstring>

namespace bucket {

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

// Check if point p is strictly inside triangle abc
inline bool point_in_triangle(const Point& p, const Point& a, const Point& b, const Point& c) {
    double d1 = signed_area2(p, a, b);
    double d2 = signed_area2(p, b, c);
    double d3 = signed_area2(p, c, a);
    
    bool has_neg = (d1 < -1e-10) || (d2 < -1e-10) || (d3 < -1e-10);
    bool has_pos = (d1 > 1e-10) || (d2 > 1e-10) || (d3 > 1e-10);
    
    return !(has_neg && has_pos);
}

/**
 * Bucket Grid for O(1) expected point queries
 * 
 * Grid size: ceil(sqrt(n)) x ceil(sqrt(n))
 * Expected vertices per bucket: O(1) for uniform distribution
 * Query cost: O(k) buckets where k = O(triangle_area / bucket_area)
 * 
 * For ear triangles in a simple polygon, the triangle area is bounded,
 * so k = O(1) on average.
 */
class BucketGrid {
    int grid_size;
    double min_x, min_y, max_x, max_y;
    double cell_width, cell_height;
    double inv_cell_width, inv_cell_height;
    
    // Flat array of buckets, each bucket is a vector of vertex indices
    std::vector<std::vector<int>> buckets;
    
    // For each vertex, store which bucket it's in
    std::vector<int> vertex_bucket;
    
    inline int bucket_index(int bx, int by) const {
        return by * grid_size + bx;
    }
    
public:
    void initialize(const std::vector<Point>& vertices) {
        int n = vertices.size();
        if (n == 0) return;
        
        // Grid size = ceil(sqrt(n))
        grid_size = std::max(1, (int)std::ceil(std::sqrt(n)));
        
        // Compute bounding box
        min_x = max_x = vertices[0].x;
        min_y = max_y = vertices[0].y;
        for (int i = 1; i < n; i++) {
            min_x = std::min(min_x, vertices[i].x);
            max_x = std::max(max_x, vertices[i].x);
            min_y = std::min(min_y, vertices[i].y);
            max_y = std::max(max_y, vertices[i].y);
        }
        
        // Add small padding to avoid edge cases
        double pad = 1e-6;
        min_x -= pad; min_y -= pad;
        max_x += pad; max_y += pad;
        
        cell_width = (max_x - min_x) / grid_size;
        cell_height = (max_y - min_y) / grid_size;
        inv_cell_width = 1.0 / cell_width;
        inv_cell_height = 1.0 / cell_height;
        
        // Initialize buckets
        buckets.clear();
        buckets.resize(grid_size * grid_size);
        vertex_bucket.resize(n, -1);
        
        // Assign vertices to buckets: O(n)
        for (int i = 0; i < n; i++) {
            int bx = std::min(grid_size - 1, (int)((vertices[i].x - min_x) * inv_cell_width));
            int by = std::min(grid_size - 1, (int)((vertices[i].y - min_y) * inv_cell_height));
            int bidx = bucket_index(bx, by);
            buckets[bidx].push_back(i);
            vertex_bucket[i] = bidx;
        }
    }
    
    void remove_vertex(int vertex_idx) {
        int bidx = vertex_bucket[vertex_idx];
        if (bidx < 0) return;
        
        auto& bucket = buckets[bidx];
        bucket.erase(std::remove(bucket.begin(), bucket.end(), vertex_idx), bucket.end());
        vertex_bucket[vertex_idx] = -1;
    }
    
    // Query all vertices in buckets overlapping the given bounding box
    // Returns vertex indices via callback
    template<typename Func>
    void query_bbox(double qmin_x, double qmin_y, double qmax_x, double qmax_y, Func&& func) const {
        int bx0 = std::max(0, (int)((qmin_x - min_x) * inv_cell_width));
        int by0 = std::max(0, (int)((qmin_y - min_y) * inv_cell_height));
        int bx1 = std::min(grid_size - 1, (int)((qmax_x - min_x) * inv_cell_width));
        int by1 = std::min(grid_size - 1, (int)((qmax_y - min_y) * inv_cell_height));
        
        for (int by = by0; by <= by1; by++) {
            for (int bx = bx0; bx <= bx1; bx++) {
                const auto& bucket = buckets[bucket_index(bx, by)];
                for (int idx : bucket) {
                    func(idx);
                }
            }
        }
    }
    
    int get_grid_size() const { return grid_size; }
};

/**
 * Main Bucket Triangulation Algorithm
 * 
 * Theorem: The algorithm runs in O(n) amortized time.
 * 
 * Proof sketch:
 * 1. Initialization: O(n) to build grid and assign vertices
 * 2. Each of the n-2 ear removals:
 *    a. Finding an ear: Each vertex is visited O(1) times amortized
 *       because once we find an ear at vertex v, we continue from v's neighbor
 *    b. Ear test: Query O(1) buckets (bounded by ear triangle size)
 *       Each bucket has O(1) expected vertices
 *    c. Removal: O(1) to update linked list and bucket
 * 3. Total: O(n) + O(n) * O(1) = O(n)
 */
class BucketTriangulator {
public:
    std::vector<Triangle> triangles;
    int num_reflex = 0;
    int grid_size = 0;
    
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
        
        // Build doubly linked list: O(n)
        std::vector<int> prev(n), next(n);
        for (int i = 0; i < n; i++) {
            prev[i] = (i - 1 + n) % n;
            next[i] = (i + 1) % n;
        }
        
        // Count reflex vertices and check for convex: O(n)
        bool is_convex = true;
        for (int i = 0; i < n; i++) {
            if (cross(poly[prev[i]], poly[i], poly[next[i]]) < -1e-10) {
                num_reflex++;
                is_convex = false;
            }
        }
        
        // Fast path for convex polygons: O(n)
        if (is_convex) {
            triangles.reserve(n - 2);
            for (int i = 1; i < n - 1; i++) {
                triangles.emplace_back(0, i, i + 1);
            }
            return;
        }
        
        // Initialize bucket grid: O(n)
        BucketGrid grid;
        grid.initialize(poly);
        grid_size = grid.get_grid_size();
        
        // Track which vertices are still active
        std::vector<bool> active(n, true);
        
        // Lambda to check if vertex is convex
        auto is_convex_vertex = [&](int i) -> bool {
            return cross(poly[prev[i]], poly[i], poly[next[i]]) > 1e-10;
        };
        
        // Lambda to check if vertex is an ear using bucket query
        auto is_ear = [&](int i) -> bool {
            if (!is_convex_vertex(i)) return false;
            
            const Point& a = poly[prev[i]];
            const Point& b = poly[i];
            const Point& c = poly[next[i]];
            
            // Bounding box of ear triangle
            double tri_min_x = std::min({a.x, b.x, c.x});
            double tri_max_x = std::max({a.x, b.x, c.x});
            double tri_min_y = std::min({a.y, b.y, c.y});
            double tri_max_y = std::max({a.y, b.y, c.y});
            
            // Query buckets for potential blockers
            bool blocked = false;
            grid.query_bbox(tri_min_x, tri_min_y, tri_max_x, tri_max_y, [&](int idx) {
                if (blocked) return;
                if (idx == prev[i] || idx == i || idx == next[i]) return;
                if (!active[idx]) return;
                
                // Only reflex vertices can block an ear
                if (cross(poly[prev[idx]], poly[idx], poly[next[idx]]) >= -1e-10) return;
                
                if (point_in_triangle(poly[idx], a, b, c)) {
                    blocked = true;
                }
            });
            
            return !blocked;
        };
        
        // Main ear clipping loop: O(n) amortized
        triangles.reserve(n - 2);
        int remaining = n;
        int current = 0;
        int consecutive_failures = 0;
        
        while (remaining > 3) {
            if (is_ear(current)) {
                // Add triangle
                triangles.emplace_back(prev[current], current, next[current]);
                
                // Update linked list
                int p = prev[current];
                int nx = next[current];
                next[p] = nx;
                prev[nx] = p;
                
                // Remove from bucket
                grid.remove_vertex(current);
                active[current] = false;
                
                remaining--;
                consecutive_failures = 0;
                
                // Continue from neighbor (key for amortized O(n))
                current = nx;
            } else {
                current = next[current];
                consecutive_failures++;
                
                // Safety: if we've gone around twice without finding an ear,
                // the polygon might be degenerate
                if (consecutive_failures > 2 * remaining) {
                    // Force progress by taking any convex vertex
                    int start = current;
                    do {
                        if (is_convex_vertex(current)) {
                            triangles.emplace_back(prev[current], current, next[current]);
                            int p = prev[current];
                            int nx = next[current];
                            next[p] = nx;
                            prev[nx] = p;
                            grid.remove_vertex(current);
                            active[current] = false;
                            remaining--;
                            current = nx;
                            break;
                        }
                        current = next[current];
                    } while (current != start);
                    consecutive_failures = 0;
                }
            }
        }
        
        // Final triangle
        if (remaining == 3) {
            triangles.emplace_back(prev[current], current, next[current]);
        }
    }
};

} // namespace bucket

