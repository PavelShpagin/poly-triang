/**
 * CLI for Reflex Triangulation Algorithm
 * 
 * Usage: reflex_cli --input polygon.poly --output triangulation.tri
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <cstring>

#include "reflex_triangulation.hpp"

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
    
    std::vector<reflex::Point> polygon(n);
    for (int i = 0; i < n; i++) {
        fin >> polygon[i].x >> polygon[i].y;
        polygon[i].index = i;
    }
    fin.close();
    
    // Triangulate with timing
    reflex::ReflexTriangulator triangulator;
    
    auto start = std::chrono::high_resolution_clock::now();
    triangulator.triangulate(polygon);
    auto end = std::chrono::high_resolution_clock::now();
    
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
    fout << triangulator.triangles.size() << "\n";
    for (const auto& tri : triangulator.triangles) {
        fout << tri.v0 << " " << tri.v1 << " " << tri.v2 << "\n";
    }
    fout.close();
    
    // Output in benchmark format
    std::cout << "reflex,vertices=" << n 
              << ",triangles=" << triangulator.triangles.size()
              << ",reflex_count=" << triangulator.num_reflex
              << ",time_ms=" << elapsed_ms << "\n";
    
    return 0;
}

