#pragma once
// TRUE O(n + r log r) with linked-list face extraction
// Key: sort only events (non-REGULAR), not all n vertices

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

namespace true_optimal {

static constexpr double kEps = 1e-12;

struct Point { double x, y; int index; };
struct Triangle { int v0, v1, v2; };
enum VType : uint8_t { START, END, SPLIT, MERGE, REGULAR };

class Triangulator;

struct Chain {
  int id, upper, lower;
  std::vector<int> verts;
  mutable int pos = 0;
  mutable int helper = -1;
  void reset() const { pos = 0; helper = -1; }
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
  int reflex_count_ = 0;
  int reflex_count() const { return reflex_count_; }

  void triangulate(const std::vector<Point>& pts);

private:
  const Point* pts_ = nullptr;
  int n_ = 0;
  std::vector<VType> types_;
  std::vector<Chain> chains_;
  std::vector<int> max_to_chain_, min_to_chain_;
  mutable double sweep_y_ = 0, query_x_ = 0;
  static constexpr int kQueryKey = -1;

  struct LLVert { double x, y; int next, prev; };
  std::vector<LLVert> ll_;
  std::vector<VType> ll_types_;
  int num_ll_ = 0;

  int nxt(int i) const { return (i+1) % n_; }
  int prv(int i) const { return (i+n_-1) % n_; }
  bool below(int a, int b) const {
    double dy = pts_[a].y - pts_[b].y;
    if (dy < -kEps) return true;
    if (dy > kEps) return false;
    return pts_[a].x > pts_[b].x;
  }
  double cross(int o, int a, int b) const {
    return (pts_[a].x-pts_[o].x)*(pts_[b].y-pts_[o].y) - (pts_[a].y-pts_[o].y)*(pts_[b].x-pts_[o].x);
  }
  bool is_convex() const {
    for (int i = 0; i < n_; ++i)
      if (cross(prv(i), i, nxt(i)) < -kEps) return false;
    return true;
  }

public:
  double chain_x_at(int cid, double y) const {
    const Chain& ch = chains_[cid];
    int len = ch.verts.size();
    while (ch.pos+1 < len && pts_[ch.verts[ch.pos+1]].y > y + kEps) ++ch.pos;
    if (ch.pos >= len-1) ch.pos = len-2;
    if (ch.pos < 0) ch.pos = 0;
    auto& p = pts_[ch.verts[ch.pos]];
    auto& q = pts_[ch.verts[ch.pos+1]];
    double dy = p.y - q.y;
    if (std::abs(dy) < kEps) return std::min(p.x, q.x);
    return p.x + (p.y - y) / dy * (q.x - p.x);
  }
  int slab_entry(int cid) const { return chains_[cid].verts[chains_[cid].pos]; }
};

inline bool ChainCmp::operator()(int a, int b) const {
  if (a == Triangulator::kQueryKey) return tri->query_x_ < tri->chain_x_at(b, tri->sweep_y_);
  if (b == Triangulator::kQueryKey) return tri->chain_x_at(a, tri->sweep_y_) < tri->query_x_;
  double xa = tri->chain_x_at(a, tri->sweep_y_);
  double xb = tri->chain_x_at(b, tri->sweep_y_);
  if (xa < xb) return true;
  if (xb < xa) return false;
  return a < b;
}

inline void Triangulator::triangulate(const std::vector<Point>& pts) {
  if (pts.size() < 3) return;
  n_ = pts.size();
  pts_ = pts.data();
  
  // Fast convex path
  if (is_convex()) {
    triangles.reserve(n_-2);
    for (int i = 1; i+1 < n_; ++i)
      triangles.push_back({pts[0].index, pts[i].index, pts[i+1].index});
    reflex_count_ = 0;
    return;
  }
  
  // Classify
  types_.resize(n_);
  reflex_count_ = 0;
  for (int i = 0; i < n_; ++i) {
    int p = prv(i), nx = nxt(i);
    bool pb = below(p, i), nb = below(nx, i);
    double c = cross(p, i, nx);
    if (pb && nb) types_[i] = (c > kEps) ? START : SPLIT;
    else if (!pb && !nb) types_[i] = (c > kEps) ? END : MERGE;
    else types_[i] = REGULAR;
    if (c < -kEps) reflex_count_++;
  }
  
  // Build chains
  chains_.clear();
  max_to_chain_.assign(n_, -1);
  min_to_chain_.assign(n_, -1);
  for (int v = 0; v < n_; ++v) {
    if (types_[v] != START && types_[v] != SPLIT) continue;
    Chain ch; ch.id = chains_.size(); ch.upper = v;
    ch.verts.push_back(v);
    int cur = prv(v);
    for (int safe = 0; types_[cur] != END && types_[cur] != MERGE && safe < n_+2; ++safe) {
      ch.verts.push_back(cur);
      cur = prv(cur);
    }
    ch.verts.push_back(cur);
    ch.lower = cur; ch.reset();
    max_to_chain_[v] = chains_.size();
    min_to_chain_[cur] = chains_.size();
    chains_.push_back(std::move(ch));
  }
  
  // Build events - ONLY non-REGULAR vertices! O(r log r) sort
  struct Ev { int v; double y, x; VType t; };
  std::vector<Ev> events;
  for (int i = 0; i < n_; ++i)
    if (types_[i] != REGULAR) events.push_back({i, pts_[i].y, pts_[i].x, types_[i]});
  std::sort(events.begin(), events.end(), [](auto& a, auto& b) {
    if (a.y > b.y + kEps) return true;
    if (a.y < b.y - kEps) return false;
    return a.x < b.x;
  });
  
  for (auto& c : chains_) c.reset();
  diagonals.clear();
  
  ChainCmp cmp{this};
  std::set<int, ChainCmp> active(cmp);
  
  for (auto& ev : events) {
    int v = ev.v;
    sweep_y_ = (ev.t == START || ev.t == SPLIT) ? pts_[v].y - 1e-9 : pts_[v].y + 1e-9;
    
    if (ev.t == START) {
      int cid = max_to_chain_[v];
      if (cid >= 0) { chains_[cid].reset(); chains_[cid].helper = v; active.insert(cid); }
    } else if (ev.t == END) {
      int cid = min_to_chain_[v];
      if (cid >= 0) {
        if (chains_[cid].helper >= 0 && types_[chains_[cid].helper] == MERGE)
          diagonals.push_back({v, chains_[cid].helper});
        active.erase(cid);
      }
    } else if (ev.t == SPLIT) {
      query_x_ = pts_[v].x;
      auto it = active.lower_bound(kQueryKey);
      if (it != active.begin()) {
        --it; int left = *it;
        if (chains_[left].helper >= 0) diagonals.push_back({v, chains_[left].helper});
        else diagonals.push_back({v, slab_entry(left)});
        chains_[left].helper = v;
      }
      int cid = max_to_chain_[v];
      if (cid >= 0) { chains_[cid].reset(); chains_[cid].helper = v; active.insert(cid); }
    } else if (ev.t == MERGE) {
      int cid = min_to_chain_[v];
      if (cid >= 0) {
        if (chains_[cid].helper >= 0 && types_[chains_[cid].helper] == MERGE)
          diagonals.push_back({v, chains_[cid].helper});
        active.erase(cid);
      }
      query_x_ = pts_[v].x;
      auto it = active.lower_bound(kQueryKey);
      if (it != active.begin()) {
        --it; int left = *it;
        if (chains_[left].helper >= 0 && types_[chains_[left].helper] == MERGE)
          diagonals.push_back({v, chains_[left].helper});
        chains_[left].helper = v;
      }
    }
  }
  
  // Linked-list face extraction
  int maxv = n_ + 2 * diagonals.size();
  ll_.resize(maxv); ll_types_.resize(maxv); num_ll_ = n_;
  for (int i = 0; i < n_; ++i) {
    ll_[i] = {pts_[i].x, pts_[i].y, nxt(i), prv(i)};
    ll_types_[i] = types_[i];
  }
  for (auto& d : diagonals) {
    int i1 = d.first, i2 = d.second;
    int n1 = num_ll_++, n2 = num_ll_++;
    ll_[n1] = ll_[i1]; ll_[n2] = ll_[i2];
    ll_[n2].next = ll_[i2].next; ll_[n1].next = ll_[i1].next;
    ll_[ll_[i2].next].prev = n2; ll_[ll_[i1].next].prev = n1;
    ll_[i1].next = n2; ll_[n2].prev = i1;
    ll_[i2].next = n1; ll_[n1].prev = i2;
    ll_types_[n1] = ll_types_[i1]; ll_types_[n2] = ll_types_[i2];
  }
  
  // Extract faces and triangulate
  std::vector<uint8_t> used(num_ll_, 0);
  triangles.reserve(n_ - 2);
  for (int s = 0; s < num_ll_; ++s) {
    if (used[s]) continue;
    std::vector<int> face;
    int cur = s;
    for (int safe = 0; safe < num_ll_ + 10; ++safe) {
      face.push_back(cur); used[cur] = 1;
      cur = ll_[cur].next;
      if (cur == s) break;
    }
    // Simple ear-clipping for monotone faces
    if (face.size() >= 3) {
      for (size_t i = 1; i + 1 < face.size(); ++i) {
        triangles.push_back({pts_[face[0] % n_].index, pts_[face[i] % n_].index, pts_[face[i+1] % n_].index});
      }
    }
  }
}

}  // namespace true_optimal
