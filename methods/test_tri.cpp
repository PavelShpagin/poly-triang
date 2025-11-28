#include <iostream>
#include <fstream>
#include <vector>
#include "include/reflex_chord.hpp"

int main() {
    // Simple test: square
    std::vector<reflex_chord::Point> square = {
        {0, 0, 0}, {1, 0, 1}, {1, 1, 2}, {0, 1, 3}
    };
    
    reflex_chord::ReflexChordTriangulator tri;
    auto result = tri.triangulate(square);
    
    std::cout << "Square (4 vertices): " << result.size() << " triangles (expected 2)" << std::endl;
    for (auto& t : result) {
        std::cout << "  (" << t.v0 << ", " << t.v1 << ", " << t.v2 << ")" << std::endl;
    }
    
    // L-shape (6 vertices, 1 reflex)
    std::vector<reflex_chord::Point> lshape = {
        {0, 0, 0}, {2, 0, 1}, {2, 1, 2}, {1, 1, 3}, {1, 2, 4}, {0, 2, 5}
    };
    result = tri.triangulate(lshape);
    std::cout << "L-shape (6 vertices): " << result.size() << " triangles (expected 4)" << std::endl;
    for (auto& t : result) {
        std::cout << "  (" << t.v0 << ", " << t.v1 << ", " << t.v2 << ")" << std::endl;
    }
    
    return 0;
}

