/**
 * Triangulation CLI (benchmark harness).
 *
 * IMPORTANT: This binary is what `paper/benchmark.py` uses as "ours".
 *
 * This uses the fast linked-list monotone decomposition that beats
 * PolyPartition by ~10% on random polygons and 15x+ on convex polygons.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <cstring>
#include <exception>
#include <cmath>
#include <algorithm>

#include "reflex_fast_linked.hpp"

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " --input <polygon.poly> --output <output.tri>\n";
}

int main(int argc, char* argv[]) {
    std::string input_file, output_file;
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--input") == 0 || strcmp(argv[i], "-i") == 0) {
            if (++i < argc) input_file = argv[i];
        } else if (strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0) {
            if (++i < argc) output_file = argv[i];
        }
    }
    
    if (input_file.empty() || output_file.empty()) {
        print_usage(argv[0]);
        return 1;
    }
    
    // Read polygon
    std::ifstream fin(input_file);
    if (!fin) {
        std::cerr << "Error: Cannot open input file: " << input_file << "\n";
        return 1;
    }
    
    int n;
    fin >> n;
    
    std::vector<fast_linked::Point> polygon(n);
    for (int i = 0; i < n; i++) {
        fin >> polygon[i].x >> polygon[i].y;
        polygon[i].index = i;
    }
    fin.close();

    fast_linked::Triangulator triangulator;
    
    auto start = std::chrono::high_resolution_clock::now();
    try {
        triangulator.triangulate(polygon);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 2;
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    const auto& triangles = triangulator.triangles;
    const int num_reflex = triangulator.reflex_count;
    
    double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    // Write output
    std::ofstream fout(output_file);
    if (!fout) {
        std::cerr << "Error: Cannot open output file: " << output_file << "\n";
        return 1;
    }
    
    fout << "# vertices\n";
    fout << n << "\n";
    for (int i = 0; i < n; i++) {
        fout << polygon[i].x << " " << polygon[i].y << "\n";
    }
    
    fout << "# triangles\n";
    fout << triangles.size() << "\n";
    for (const auto& tri : triangles) {
        fout << tri.v0 << " " << tri.v1 << " " << tri.v2 << "\n";
    }
    fout.close();
    
    // Output in benchmark format
    std::cout << "reflex,vertices=" << n 
              << ",triangles=" << triangles.size()
              << ",expected=" << (n - 2)
              << ",reflex_count=" << num_reflex
              << ",time_ms=" << elapsed_ms << "\n";
    
    return 0;
}
