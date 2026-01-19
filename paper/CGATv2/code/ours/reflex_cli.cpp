/**
 * Triangulation CLI (CGAT artifact / benchmark harness).
 *
 * Used by `paper/CGAT/tools/benchmark_cgat.py`.
 *
 * Algorithms:
 * - chain_only: chain-based sweep (no heuristics / no fallback).
 * - chain: alias of chain_only (kept for backwards compatibility).
 * - linked: edge-based sweep with a linked representation.
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
#include <time.h>

#include "reflex_fast_linked.hpp"
#include "reflex_chain_triangulate.hpp"

// k = number of local maxima (equivalently, local minima) w.r.t. the sweep direction.
// This is the event-complexity parameter used in the paper (O(n + k log k)).
static inline bool below_k(const std::vector<std::pair<double, double>>& coords, int a, int b) {
    const double ay = coords[a].second;
    const double by = coords[b].second;
    if (ay < by) return true;
    if (ay > by) return false;
    const double ax = coords[a].first;
    const double bx = coords[b].first;
    if (ax > bx) return true;
    if (ax < bx) return false;
    return a > b;
}

static int count_local_maxima_k(const std::vector<std::pair<double, double>>& coords) {
    const int n = static_cast<int>(coords.size());
    if (n < 3) return 0;
    int k = 0;
    for (int i = 0; i < n; ++i) {
        const int p = (i - 1 + n) % n;
        const int nx = (i + 1) % n;
        if (below_k(coords, p, i) && below_k(coords, nx, i)) ++k;
    }
    return k;
}

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " --input <polygon.poly> --output <output.tri> [--algo chain|chain_only|linked]\n";
}

int main(int argc, char* argv[]) {
    std::string input_file, output_file;
    std::string algo = "chain_only"; // core method (no heuristics / no fallback)
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--input") == 0 || strcmp(argv[i], "-i") == 0) {
            if (++i < argc) input_file = argv[i];
        } else if (strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0) {
            if (++i < argc) output_file = argv[i];
        } else if (strcmp(argv[i], "--algo") == 0 || strcmp(argv[i], "-a") == 0) {
            if (++i < argc) algo = argv[i];
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
    
    std::vector<std::pair<double, double>> coords(n);
    for (int i = 0; i < n; i++) {
        fin >> coords[i].first >> coords[i].second;
    }
    fin.close();
    const int k_count = count_local_maxima_k(coords);

    if (algo == "linked") {
        std::vector<fast_linked::Point> polygon(n);
        for (int i = 0; i < n; i++) {
            polygon[i].x = coords[i].first;
            polygon[i].y = coords[i].second;
            polygon[i].index = i;
        }

        fast_linked::Triangulator triangulator;
        auto cpu_now_ms = []() -> double {
            struct timespec ts;
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
            return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
        };
        const double t0 = cpu_now_ms();
        try {
            triangulator.triangulate(polygon);
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << "\n";
            return 2;
        }
        const double t1 = cpu_now_ms();

        const auto& triangles = triangulator.triangles;
        const int num_reflex = triangulator.reflex_count;
        double elapsed_ms = (t1 - t0);

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

        std::cout << "reflex,mode=linked,vertices=" << n
                  << ",triangles=" << triangles.size()
                  << ",expected=" << (n - 2)
                  << ",reflex_count=" << num_reflex
                  << ",k_count=" << k_count
                  << ",time_ms=" << elapsed_ms << "\n";
        return 0;
    }

    if (algo == "chain" || algo == "chain_only") {
        std::vector<reflex_tri::Point> polygon(n);
        for (int i = 0; i < n; i++) {
            polygon[i].x = coords[i].first;
            polygon[i].y = coords[i].second;
            polygon[i].index = i;
        }

        reflex_tri::Triangulator triangulator;

        auto cpu_now_ms = []() -> double {
            struct timespec ts;
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
            return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
        };
        const double t0 = cpu_now_ms();
        try {
            triangulator.triangulate(polygon);
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << "\n";
            return 2;
        }
        const double t1 = cpu_now_ms();

        const auto& triangles = triangulator.debug_triangles();
        const int num_reflex = triangulator.reflex_count();
        double elapsed_ms = (t1 - t0);

        // NOTE: chain triangulator may reverse vertices internally to ensure CCW;
        // we therefore write the (possibly reordered) vertices we triangulated.
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

        std::cout << "reflex,mode=chain_only,vertices=" << n
                  << ",triangles=" << triangles.size()
                  << ",expected=" << (n - 2)
                  << ",reflex_count=" << num_reflex
                  << ",k_count=" << k_count
                  << ",time_ms=" << elapsed_ms << "\n";
        return 0;
    }

    std::cerr << "Error: unknown --algo '" << algo << "' (expected 'chain_only', 'chain', or 'linked')\n";
    return 1;
}
