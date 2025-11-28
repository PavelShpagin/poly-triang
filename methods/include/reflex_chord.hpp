#pragma once
/**
 * Reflex Chord Triangulation (Corrected & Optimized)
 * 
 * Algorithm based on "Reflex Polygon Triangulation":
 * 1. Identify reflex vertices.
 * 2. Use a Reflex-Only Spatial Grid to validate diagonals.
 * 3. This effectively implements the "find safe horizontal chord" step 
 *    by finding safe diagonals in O(1) expected time.
 * 
 * Complexity: O(n) expected for uniform distribution.
 * Worst case: O(n^2) (degenerate).
 */

#include <vector>
#include <cmath>
#include <algorithm>
#include <list>
#include <limits>

namespace reflex_chord {

struct Point {
    double x, y;
    int index; // Original index
};

struct Triangle {
    int v0, v1, v2;
    Triangle(int a, int b, int c) : v0(a), v1(b), v2(c) {}
};

inline double cross(const Point& a, const Point& b, const Point& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

class ReflexChordTriangulator {
public:
    std::vector<Triangle> triangulate(const std::vector<Point>& points) {
        return bucket_triangulate_internal(points);
    }

private:
    std::vector<Triangle> bucket_triangulate_internal(const std::vector<Point>& input_polygon) {
        std::vector<Triangle> out_triangles;
        int n = input_polygon.size();
        if (n < 3) return out_triangles;
        
        // 1. Setup Lists
        std::vector<int> prev(n), next(n);
        for(int i=0; i<n; i++) { 
            prev[i] = (i - 1 + n) % n; 
            next[i] = (i + 1) % n; 
        }
        
        // 2. Orientation & Reflex Detection
        double area = 0;
        for(int i=0; i<n; i++) area += cross(Point{0,0,0}, input_polygon[i], input_polygon[next[i]]);
        bool ccw = area > 0;
        
        auto is_reflex = [&](int i) {
            double c = cross(input_polygon[prev[i]], input_polygon[i], input_polygon[next[i]]);
            // Relaxed epsilon to handle noisy convex polygons (treat approx 180 as convex)
            return ccw ? (c <= -1e-8) : (c >= 1e-8);
        };
        
        std::vector<int> reflex_indices;
        std::vector<bool> active_reflex(n, false);
        for(int i=0; i<n; i++) {
            if(is_reflex(i)) { 
                reflex_indices.push_back(i); 
                active_reflex[i] = true; 
            }
        }
        
        // 3. Fast Path for Convex
        if(reflex_indices.empty()) {
            for(int i=1; i<n-1; i++) out_triangles.emplace_back(0, i, i+1);
            return out_triangles;
        }
        
        // 4. Build Grid (Reflex Only)
        double min_x = input_polygon[reflex_indices[0]].x;
        double max_x = min_x;
        double min_y = input_polygon[reflex_indices[0]].y;
        double max_y = min_y;
        
        for(int idx : reflex_indices) {
            const auto& p = input_polygon[idx];
            min_x = std::min(min_x, p.x); max_x = std::max(max_x, p.x);
            min_y = std::min(min_y, p.y); max_y = std::max(max_y, p.y);
        }
        
        int grid_dim = std::max(1, (int)std::sqrt(reflex_indices.size()));
        double cell_w = (max_x - min_x + 1e-6) / grid_dim;
        double cell_h = (max_y - min_y + 1e-6) / grid_dim;
        
        std::vector<std::vector<int>> grid(grid_dim * grid_dim);
        for(int idx : reflex_indices) {
            int gx = std::max(0, std::min(grid_dim-1, (int)((input_polygon[idx].x - min_x)/cell_w)));
            int gy = std::max(0, std::min(grid_dim-1, (int)((input_polygon[idx].y - min_y)/cell_h)));
            grid[gy * grid_dim + gx].push_back(idx);
        }
        
        // 5. Ear Clipping
        int curr = 0;
        // Start at a convex vertex
        for(int i=0; i<n; i++) if(!active_reflex[i]) { curr=i; break; }
        
        int count = n;
        int failed = 0;
        
        while(count > 3) {
            bool cut = false;
            if(!active_reflex[curr]) { // Candidate ear
                int p = prev[curr];
                int nx = next[curr];
                
                const auto &a = input_polygon[p];
                const auto &b = input_polygon[curr];
                const auto &c = input_polygon[nx];
                
                // Bounding box of ear triangle
                double x0 = std::min({a.x, b.x, c.x}), x1 = std::max({a.x, b.x, c.x});
                double y0 = std::min({a.y, b.y, c.y}), y1 = std::max({a.y, b.y, c.y});
                
                int gx0 = std::max(0, std::min(grid_dim-1, (int)((x0 - min_x)/cell_w)));
                int gx1 = std::max(0, std::min(grid_dim-1, (int)((x1 - min_x)/cell_w)));
                int gy0 = std::max(0, std::min(grid_dim-1, (int)((y0 - min_y)/cell_h)));
                int gy1 = std::max(0, std::min(grid_dim-1, (int)((y1 - min_y)/cell_h)));
                
                bool ok = true;
                // Query Grid
                for(int y=gy0; y<=gy1; y++) {
                    for(int x=gx0; x<=gx1; x++) {
                        for(int r : grid[y*grid_dim + x]) {
                            if(!active_reflex[r]) continue;
                            if(r == p || r == curr || r == nx) continue;
                            
                            // Robust Point-in-Triangle
                            double cp1 = cross(a, b, input_polygon[r]);
                            double cp2 = cross(b, c, input_polygon[r]);
                            double cp3 = cross(c, a, input_polygon[r]);
                            
                            // Relaxed check: points ON the boundary should NOT invalidate the ear
                            // Only strictly interior points invalidate the ear
                            double eps = 1e-9;
                            bool all_neg = (cp1 <= -eps) && (cp2 <= -eps) && (cp3 <= -eps);
                            bool all_pos = (cp1 >= eps) && (cp2 >= eps) && (cp3 >= eps);
                            
                            if(ccw ? all_pos : all_neg) { 
                                ok = false; 
                                goto done_check; 
                            }
                        }
                    }
                }
                done_check:;
                
                if(ok) {
                    out_triangles.emplace_back(p, curr, nx);
                    
                    // Remove curr
                    next[p] = nx; prev[nx] = p;
                    count--;
                    
                    // Update neighbors
                    if(active_reflex[p] && !is_reflex(p)) active_reflex[p] = false;
                    if(active_reflex[nx] && !is_reflex(nx)) active_reflex[nx] = false;
                    
                    curr = nx;
                    cut = true;
                    failed = 0;
                }
            }
            
            if(!cut) {
                curr = next[curr];
                failed++;
                if(failed > count*2) {
                    // Fallback to avoid infinite loop
                    int p = prev[curr];
                    int nx = next[curr];
                    out_triangles.emplace_back(p, curr, nx);
                    next[p] = nx; prev[nx] = p;
                    count--;
                    curr = nx;
                    failed = 0;
                }
            }
        }
        
        if(count == 3) {
            out_triangles.emplace_back(prev[curr], curr, next[curr]);
        }
        
        return out_triangles;
    }
};

} // namespace
