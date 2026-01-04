#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <vector>

// Copy of the triangulate_monotone to debug
static constexpr double kEps = 1e-12;

struct Vertex {
  double x, y;
  int next, prev;
};

int triangulate_monotone_count(const std::vector<int>& face, const std::vector<Vertex>& verts, int n) {
  const int m = static_cast<int>(face.size());
  
  if (m < 3) return 0;
  if (m == 3) return 1;
  
  // Find top vertex
  int top = 0;
  for (int i = 1; i < m; ++i) {
    if (verts[face[i]].y > verts[face[top]].y + kEps ||
        (std::abs(verts[face[i]].y - verts[face[top]].y) < kEps &&
         verts[face[i]].x < verts[face[top]].x)) {
      top = i;
    }
  }
  
  // Sort by y descending
  std::vector<std::pair<int, int>> sorted(m);
  for (int i = 0; i < m; ++i) {
    sorted[i] = {i, face[i]};
  }
  std::sort(sorted.begin(), sorted.end(), [&verts](const auto& a, const auto& b) {
    double ya = verts[a.second].y, yb = verts[b.second].y;
    if (std::abs(ya - yb) > kEps) return ya > yb;
    return verts[a.second].x < verts[b.second].x;
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
  
  int tri_count = 0;
  std::vector<int> stk;
  stk.push_back(sorted[0].first);
  stk.push_back(sorted[1].first);
  
  for (int i = 2; i < m; ++i) {
    int vi = sorted[i].first;
    
    if (side[vi] != side[stk.back()]) {
      while (stk.size() > 1) {
        tri_count++;
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
        double c = (verts[face[a]].x - verts[face[vi]].x) *
                   (verts[face[last]].y - verts[face[vi]].y) -
                   (verts[face[a]].y - verts[face[vi]].y) *
                   (verts[face[last]].x - verts[face[vi]].x);
        bool ok = (side[vi] == 0) ? (c > kEps) : (c < -kEps);
        if (!ok) break;
        tri_count++;
        last = stk.back();
        stk.pop_back();
      }
      stk.push_back(last);
      stk.push_back(vi);
    }
  }
  
  return tri_count;
}

int main() {
  // Test with a simple monotone polygon
  std::vector<Vertex> verts = {
    {0, 4, 1, 4},  // 0: top
    {1, 3, 2, 0},  // 1
    {2, 2, 3, 1},  // 2
    {3, 1, 4, 2},  // 3
    {2, 0, 0, 3},  // 4: bottom
  };
  
  std::vector<int> face = {0, 1, 2, 3, 4};
  int expected = static_cast<int>(face.size()) - 2;
  int actual = triangulate_monotone_count(face, verts, 5);
  
  std::cout << "Simple monotone: m=" << face.size() << " expected=" << expected << " actual=" << actual << std::endl;
  
  // Test with varying sizes
  for (int m = 3; m <= 10; ++m) {
    std::vector<Vertex> v(m);
    std::vector<int> f(m);
    for (int i = 0; i < m; ++i) {
      v[i].x = i;
      v[i].y = m - i;
      v[i].next = (i + 1) % m;
      v[i].prev = (i - 1 + m) % m;
      f[i] = i;
    }
    int exp = m - 2;
    int act = triangulate_monotone_count(f, v, m);
    std::cout << "m=" << m << " expected=" << exp << " actual=" << act;
    if (exp != act) std::cout << " MISMATCH!";
    std::cout << std::endl;
  }
  
  return 0;
}

