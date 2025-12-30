#include <iostream>
#include <vector>
#include "reflex_optimal.hpp"

int main() {
    // Paper example polygon (CCW)
    std::vector<reflex_opt::Point> poly = {
        {0.0, 2.5, 0},
        {1.0, 5.5, 1},
        {2.5, 4.0, 2},
        {4.0, 6.5, 3},
        {5.5, 5.0, 4},
        {7.0, 7.0, 5},
        {8.0, 5.5, 6},
        {6.5, 3.5, 7},
        {8.0, 1.5, 8},
        {5.0, 2.5, 9},
        {3.0, 0.0, 10},
        {1.5, 1.5, 11},
    };
    
    reflex_opt::Triangulator tri;
    const auto tris = tri.triangulate(poly);
    
    std::cout << "Paper example (12 vertices):" << std::endl;
    std::cout << "  Reflex vertices r: " << tri.reflex_count() << std::endl;
    std::cout << "  Triangles: " << tris.size() << " (expected 10)" << std::endl;
    
    return 0;
}

