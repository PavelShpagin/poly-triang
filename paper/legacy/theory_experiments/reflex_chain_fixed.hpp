#pragma once
/**
 * TRUE O(n + r log r) Polygon Triangulation - FIXED
 * 
 * Key fixes:
 * 1. Chains go via nxt() not prv() - LEFT boundary when traversing CCW
 * 2. Proper handling of both left and right chains at split/merge vertices
 * 3. Correct helper management
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

namespace reflex_chain_fixed {

static constexpr double kEps = 1e-12;

struct Point { double x, y; int index; };
struct Triangle { int v0, v1, v2; };

enum VType : uint8_t { START, END, SPLIT, MERGE, REGULAR };

class Triangulator;

struct Chain {
  int id;
  int upper, lower;        // local max -> local min
  std::vector<int> verts;  // vertices from upper to lower
  mutable int pos = 0;     // for lazy advancement
  int helper = -1;         // helper vertex for diagonal insertion
  
  void reset() { pos = 0; helper = upper; }
};

struct ChainCmp {
  const Triangulator* tri;
  bool operator()(int a, int b) const;
};

class Triangulator {
  friend struct ChainCmp;
  
public:
  std::vector<Triangle> triangles;
  std::vector<std::pair<int,int>> diagonals;
  int reflex_count = 0;
  
  void triangulate(const std::vector<Point>& input) {
    triangles.clear();
    diagonals.clear();
    reflex_count = 0;
    
    n_ = static_cast<int>(input.size());
    if (n_ < 3) return;
    
    pts_ = input;
    
    // Ensure CCW orientation
    double area = 0;
    for (int i = 0; i < n_; ++i) {
      int j = (i + 1) % n_;
      area += pts_[i].x * pts_[j].y - pts_[j].x * pts_[i].y;
    }
    if (area < 0) {
      std::reverse(pts_.begin(), pts_.end());
    }
    for (int i = 0; i < n_; ++i) pts_[i].index = i;
    
    if (n_ == 3) {
      triangles.push_back({0, 1, 2});
      return;
    }
    
    classify();
    
    // Convex fast-path
    if (reflex_count == 0) {
      for (int i = 1; i + 1 < n_; ++i)
        triangles.push_back({0, i, i + 1});
      return;
    }
    
    build_chains();
    decompose();
    extract_and_triangulate();
  }

private:
  int n_ = 0;
  std::vector<Point> pts_;
  std::vector<VType> types_;
  std::vector<Chain> chains_;
  std::vector<int> start_chain_;  // chain starting at each START/SPLIT vertex
  std::vector<int> end_chain_;    // chain ending at each END/MERGE vertex
  mutable double sweep_y_ = 0;
  mutable double query_x_ = 0;
  static constexpr int kQueryKey = -1;
  
  int nxt(int i) const { return (i + 1) % n_; }
  int prv(int i) const { return (i + n_ - 1) % n_; }
  
  bool below(int a, int b) const {
    double dy = pts_[a].y - pts_[b].y;
    if (dy < -kEps) return true;
    if (dy > kEps) return false;
    return pts_[a].x > pts_[b].x;
  }
  
  double cross(int o, int a, int b) const {
    return (pts_[a].x - pts_[o].x) * (pts_[b].y - pts_[o].y)
         - (pts_[a].y - pts_[o].y) * (pts_[b].x - pts_[o].x);
  }
  
  void classify() {
    types_.resize(n_);
    reflex_count = 0;
    
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
      
      if (c < -kEps) reflex_count++;
    }
  }
  
  void build_chains() {
    chains_.clear();
    start_chain_.assign(n_, -1);
    end_chain_.assign(n_, -1);
    
    // For each local maximum (START or SPLIT), build a chain going DOWN via nxt()
    // This chain represents the LEFT boundary of the polygon interior
    for (int v = 0; v < n_; ++v) {
      if (types_[v] != START && types_[v] != SPLIT) continue;
      
      Chain ch;
      ch.id = static_cast<int>(chains_.size());
      ch.upper = v;
      ch.verts.push_back(v);
      
      int cur = nxt(v);
      int safety = 0;
      while (types_[cur] != END && types_[cur] != MERGE && safety++ < n_ + 2) {
        ch.verts.push_back(cur);
        cur = nxt(cur);
      }
      ch.verts.push_back(cur);
      ch.lower = cur;
      ch.reset();
      
      start_chain_[v] = ch.id;
      end_chain_[cur] = ch.id;
      chains_.push_back(std::move(ch));
    }
  }
  
public:
  double chain_x_at(int cid, double y) const {
    const Chain& ch = chains_[cid];
    int len = static_cast<int>(ch.verts.size());
    
    // Lazy advancement
    while (ch.pos + 1 < len && pts_[ch.verts[ch.pos + 1]].y > y + kEps) {
      ++ch.pos;
    }
    if (ch.pos >= len - 1) ch.pos = len - 2;
    if (ch.pos < 0) ch.pos = 0;
    
    const Point& p = pts_[ch.verts[ch.pos]];
    const Point& q = pts_[ch.verts[ch.pos + 1]];
    double dy = p.y - q.y;
    if (std::abs(dy) < kEps) return std::min(p.x, q.x);
    return p.x + (p.y - y) / dy * (q.x - p.x);
  }
  
  int slab_vertex(int cid) const {
    const Chain& ch = chains_[cid];
    // Return the vertex at current slab position (upper vertex of current segment)
    return ch.verts[ch.pos];
  }

private:
  void decompose() {
    // Build events - ONLY local extrema
    struct Event { int v; double y, x; VType type; };
    std::vector<Event> events;
    
    for (int i = 0; i < n_; ++i) {
      if (types_[i] != REGULAR) {
        events.push_back({i, pts_[i].y, pts_[i].x, types_[i]});
      }
    }
    
    // O(r log r) sort
    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
      if (a.y > b.y + kEps) return true;
      if (a.y < b.y - kEps) return false;
      return a.x < b.x;
    });
    
    for (auto& c : chains_) c.reset();
    
    ChainCmp cmp{this};
    std::set<int, ChainCmp> active(cmp);
    
    for (const Event& ev : events) {
      int v = ev.v;
      
      // Set sweep_y_ slightly below for START/SPLIT (after insertion)
      // or slightly above for END/MERGE (before removal)
      if (ev.type == START || ev.type == SPLIT) {
        sweep_y_ = pts_[v].y - 1e-9;
      } else {
        sweep_y_ = pts_[v].y + 1e-9;
      }
      
      switch (ev.type) {
        case START: {
          // Insert the chain starting at v
          int cid = start_chain_[v];
          if (cid >= 0) {
            chains_[cid].reset();
            active.insert(cid);
          }
          break;
        }
        
        case END: {
          // Remove the chain ending at v
          int cid = end_chain_[v];
          if (cid >= 0) {
            // If helper is a MERGE vertex, add diagonal
            int h = chains_[cid].helper;
            if (h >= 0 && h != v && types_[h] == MERGE) {
              diagonals.push_back({v, h});
            }
            active.erase(cid);
          }
          break;
        }
        
        case SPLIT: {
          // Find left chain and add diagonal
          query_x_ = pts_[v].x;
          auto it = active.lower_bound(kQueryKey);
          if (it != active.begin()) {
            --it;
            int left_cid = *it;
            
            // Connect to helper or slab vertex
            int h = chains_[left_cid].helper;
            if (h >= 0) {
              diagonals.push_back({v, h});
            }
            chains_[left_cid].helper = v;
          }
          
          // Insert the new chain starting at v
          int cid = start_chain_[v];
          if (cid >= 0) {
            chains_[cid].reset();
            active.insert(cid);
          }
          break;
        }
        
        case MERGE: {
          // Handle the chain ending at v (right chain from interior's perspective)
          int cid = end_chain_[v];
          if (cid >= 0) {
            int h = chains_[cid].helper;
            if (h >= 0 && h != v && types_[h] == MERGE) {
              diagonals.push_back({v, h});
            }
            active.erase(cid);
          }
          
          // Find left chain and update helper
          query_x_ = pts_[v].x;
          auto it = active.lower_bound(kQueryKey);
          if (it != active.begin()) {
            --it;
            int left_cid = *it;
            int h = chains_[left_cid].helper;
            if (h >= 0 && h != v && types_[h] == MERGE) {
              diagonals.push_back({v, h});
            }
            chains_[left_cid].helper = v;
          }
          break;
        }
        
        default:
          break;
      }
    }
  }
  
  void extract_and_triangulate() {
    // Use linked-list face extraction
    int maxv = n_ + 2 * static_cast<int>(diagonals.size());
    
    struct LLV { double x, y; int next, prev; };
    std::vector<LLV> ll(maxv);
    int num_ll = n_;
    
    for (int i = 0; i < n_; ++i) {
      ll[i] = {pts_[i].x, pts_[i].y, nxt(i), prv(i)};
    }
    
    // Add diagonals
    for (auto& d : diagonals) {
      int i1 = d.first, i2 = d.second;
      int n1 = num_ll++, n2 = num_ll++;
      
      ll[n1] = ll[i1];
      ll[n2] = ll[i2];
      
      ll[n2].next = ll[i2].next;
      ll[n1].next = ll[i1].next;
      
      ll[ll[i2].next].prev = n2;
      ll[ll[i1].next].prev = n1;
      
      ll[i1].next = n2;
      ll[n2].prev = i1;
      
      ll[i2].next = n1;
      ll[n1].prev = i2;
    }
    
    // Extract faces
    std::vector<uint8_t> used(num_ll, 0);
    triangles.reserve(n_ - 2);
    
    for (int s = 0; s < num_ll; ++s) {
      if (used[s]) continue;
      
      std::vector<int> face;
      int cur = s;
      for (int safety = 0; safety < num_ll + 10; ++safety) {
        face.push_back(cur);
        used[cur] = 1;
        cur = ll[cur].next;
        if (cur == s) break;
      }
      
      if (face.size() >= 3) {
        // Fan triangulation for monotone faces
        for (size_t i = 1; i + 1 < face.size(); ++i) {
          int v0 = face[0] % n_;
          int v1 = face[i] % n_;
          int v2 = face[i + 1] % n_;
          triangles.push_back({pts_[v0].index, pts_[v1].index, pts_[v2].index});
        }
      }
    }
  }
};

inline bool ChainCmp::operator()(int a, int b) const {
  if (a == b) return false;
  double xa = (a == Triangulator::kQueryKey) ? tri->query_x_ : tri->chain_x_at(a, tri->sweep_y_);
  double xb = (b == Triangulator::kQueryKey) ? tri->query_x_ : tri->chain_x_at(b, tri->sweep_y_);
  if (xa < xb) return true;
  if (xb < xa) return false;
  return a < b;
}

}  // namespace reflex_chain_fixed

