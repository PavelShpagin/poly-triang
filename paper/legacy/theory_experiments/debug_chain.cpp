#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <set>
#include "reflex_chain_v2.hpp"
#include "../baselines/polypartition/src/polypartition.h"
#include "../baselines/polypartition/src/polypartition.cpp"

std::vector<reflex_chain_fixed::Point> make_random(int n, uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> dist(0.0, 1000.0);
  std::vector<double> angles(n);
  for (int i = 0; i < n; ++i) angles[i] = dist(rng);
  std::sort(angles.begin(), angles.end());
  std::vector<reflex_chain_fixed::Point> pts(n);
  for (int i = 0; i < n; ++i) {
    double a = angles[i] * 2.0 * M_PI / 1000.0;
    double r = 100.0 + dist(rng) * 0.5;
    pts[i].x = 500.0 + r * std::cos(a);
    pts[i].y = 500.0 + r * std::sin(a);
    pts[i].index = i;
  }
  return pts;
}

int main() {
  int n = 20, seed = 0;
  auto pts = make_random(n, seed);
  
  std::cout << "n=" << n << ", seed=" << seed << "\n\n";
  
  // Print vertices with classification
  std::cout << "Vertices (sorted by y desc):\n";
  std::vector<int> sorted(n);
  for (int i = 0; i < n; ++i) sorted[i] = i;
  std::sort(sorted.begin(), sorted.end(), [&](int a, int b) {
    if (pts[a].y != pts[b].y) return pts[a].y > pts[b].y;
    return pts[a].x < pts[b].x;
  });
  
  auto classify = [&](int i) -> const char* {
    int p = (i + n - 1) % n, nx = (i + 1) % n;
    auto below = [&](int a, int b) {
      if (pts[a].y < pts[b].y - 1e-9) return true;
      if (pts[a].y > pts[b].y + 1e-9) return false;
      return pts[a].x > pts[b].x;
    };
    bool pb = below(p, i), nb = below(nx, i);
    double c = (pts[i].x - pts[p].x) * (pts[nx].y - pts[p].y) 
             - (pts[i].y - pts[p].y) * (pts[nx].x - pts[p].x);
    if (pb && nb) return c > 0 ? "START" : "SPLIT";
    if (!pb && !nb) return c > 0 ? "END" : "MERGE";
    return "REG";
  };
  
  for (int i : sorted) {
    std::cout << "  " << i << ": (" << pts[i].x << ", " << pts[i].y << ") " << classify(i) << "\n";
  }
  
  // Our algorithm with debug
  std::cout << "\n--- Building chains ---\n";
  // Manually trace chain building
  for (int v = 0; v < n; ++v) {
    const char* t = classify(v);
    if (strcmp(t, "START") == 0 || strcmp(t, "SPLIT") == 0) {
      std::cout << "Chain from " << v << " (" << t << "): " << v;
      int cur = (v + 1) % n;
      for (int safety = 0; safety < n; ++safety) {
        const char* ct = classify(cur);
        std::cout << " -> " << cur;
        if (strcmp(ct, "END") == 0 || strcmp(ct, "MERGE") == 0) {
          std::cout << " (" << ct << ")\n";
          break;
        }
        cur = (cur + 1) % n;
      }
    }
  }
  
  reflex_chain_fixed::Triangulator tri;
  tri.triangulate(pts);
  
  std::cout << "\nOur diagonals (" << tri.diagonals.size() << "):\n";
  for (auto& d : tri.diagonals) {
    std::cout << "  " << d.first << " - " << d.second << "\n";
  }
  std::cout << "Our triangles: " << tri.triangles.size() << " (expected " << (n-2) << ")\n";
  
  // Print face details - manually extract faces
  std::cout << "\nOur faces:\n";
  {
    int maxv = n + 2 * static_cast<int>(tri.diagonals.size());
    struct LLV { double x, y; int next, prev, orig; };
    std::vector<LLV> ll(maxv);
    int num_ll = n;
    
    for (int i = 0; i < n; ++i) {
      ll[i] = {pts[i].x, pts[i].y, (i+1)%n, (i+n-1)%n, i};
    }
    
    for (auto& d : tri.diagonals) {
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
    
    std::vector<uint8_t> used(num_ll, 0);
    int fi = 0;
    for (int s = 0; s < num_ll; ++s) {
      if (used[s]) continue;
      
      std::vector<int> face;
      int cur = s;
      for (int safety = 0; safety < num_ll + 10; ++safety) {
        face.push_back(ll[cur].orig);
        used[cur] = 1;
        cur = ll[cur].next;
        if (cur == s) break;
      }
      
      std::cout << "  Face " << fi++ << " (" << face.size() << " verts): ";
      for (int v : face) std::cout << v << " ";
      std::cout << "\n";
    }
  }
  
  // PolyPartition - get monotone pieces first
  TPPLPoly poly;
  poly.Init(n);
  for (int i = 0; i < n; ++i) {
    poly[i].x = pts[i].x;
    poly[i].y = pts[i].y;
  }
  poly.SetOrientation(TPPL_ORIENTATION_CCW);
  
  std::list<TPPLPoly> inpolys;
  inpolys.push_back(poly);
  std::list<TPPLPoly> mono_pieces;
  TPPLPartition part;
  int ret = part.MonotonePartition(&inpolys, &mono_pieces);
  
  std::cout << "\nPolyPartition monotone pieces: " << mono_pieces.size() << "\n";
  
  // Print mono pieces
  std::cout << "Monotone pieces:\n";
  int pi = 0;
  for (auto& piece : mono_pieces) {
    int m = piece.GetNumPoints();
    std::cout << "  Piece " << pi++ << " (" << m << " verts): ";
    for (int i = 0; i < m; ++i) {
      // Find original index by position match
      double px = piece[i].x, py = piece[i].y;
      for (int j = 0; j < n; ++j) {
        if (std::abs(pts[j].x - px) < 0.01 && std::abs(pts[j].y - py) < 0.01) {
          std::cout << j << " ";
          break;
        }
      }
    }
    std::cout << "\n";
  }
  
  std::cout << "Implied diagonals: " << (mono_pieces.size() - 1) << "\n";
  
  // Triangulate
  std::list<TPPLPoly> triangles;
  ret = part.Triangulate_MONO(&poly, &triangles);
  std::cout << "PolyPartition triangles: " << triangles.size() << " (expected " << (n-2) << ")\n";
  
  return 0;
}

