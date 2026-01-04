#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include "../baselines/polypartition/src/polypartition.h"

using namespace std;

vector<pair<double,double>> make_random(int n, uint32_t seed) {
  mt19937 rng(seed);
  uniform_real_distribution<double> dist(0.0, 1000.0);
  vector<double> angles(n);
  for (int i = 0; i < n; ++i) angles[i] = dist(rng);
  sort(angles.begin(), angles.end());
  vector<pair<double,double>> pts(n);
  for (int i = 0; i < n; ++i) {
    double a = angles[i] * 2.0 * M_PI / 1000.0;
    double r = 100.0 + dist(rng) * 0.5;
    pts[i].first = 500.0 + r * cos(a);
    pts[i].second = 500.0 + r * sin(a);
  }
  return pts;
}

int main() {
  int n = 100;
  auto pts = make_random(n, 42);
  
  TPPLPoly poly;
  poly.Init(n);
  for (int i = 0; i < n; ++i) {
    poly[i].x = pts[i].first;
    poly[i].y = pts[i].second;
  }
  poly.SetOrientation(TPPL_ORIENTATION_CCW);
  
  list<TPPLPoly> polys, monotone, triangles;
  polys.push_back(poly);
  
  TPPLPartition pp;
  
  // First just get monotone pieces
  int result = pp.MonotonePartition(&polys, &monotone);
  cout << "MonotonePartition: " << monotone.size() << " pieces" << endl;
  
  // Count total vertices in monotone pieces
  int total_verts = 0;
  for (auto& m : monotone) {
    total_verts += m.GetNumPoints();
    cout << "  Piece: " << m.GetNumPoints() << " vertices" << endl;
  }
  
  // Now triangulate
  result = pp.Triangulate_MONO(&polys, &triangles);
  cout << "Triangulate_MONO: " << triangles.size() << " triangles (expected " << n-2 << ")" << endl;
  
  return 0;
}

