#include <algorithm>
#include "polygon_io.hpp"

#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "polytri/polytri.hpp"

using Clock = std::chrono::high_resolution_clock;

namespace {

struct Args {
    std::string input;
    std::string output;
};

Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        std::string token = argv[i];
        if ((token == "-i" || token == "--input") && i + 1 < argc) {
            args.input = argv[++i];
        } else if ((token == "-o" || token == "--output") && i + 1 < argc) {
            args.output = argv[++i];
        } else {
            std::ostringstream oss;
            oss << "Unknown or incomplete argument: " << token;
            throw std::runtime_error(oss.str());
        }
    }

    if (args.input.empty()) {
        throw std::runtime_error("Missing --input argument");
    }
    if (args.output.empty()) {
        throw std::runtime_error("Missing --output argument");
    }
    return args;
}

double signed_area(const baselines::Polygon& polygon) {
    double area = 0.0;
    const std::size_t n = polygon.size();
    for (std::size_t i = 0; i < n; ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % n];
        area += a.x * b.y - b.x * a.y;
    }
    return area * 0.5;
}

} // namespace

int main(int argc, char** argv) {
    using namespace baselines;
    try {
        const auto args = parse_args(argc, argv);
        auto polygon = read_polygon(args.input);

        if (signed_area(polygon) < 0.0) {
            std::reverse(polygon.begin(), polygon.end());
        }

        std::vector<std::vector<PolyTri::vertex_t>> contours(1);
        contours[0].reserve(polygon.size());
        for (const auto& p : polygon) {
            contours[0].emplace_back(p.x, p.y);
        }

        const auto start = Clock::now();
        auto indices = PolyTri::Triangulate(contours);
        const auto stop = Clock::now();
        const double elapsed_ms =
            std::chrono::duration<double, std::milli>(stop - start).count();

        Triangles tris;
        tris.reserve(indices.size() / 3);
        for (std::size_t i = 0; i + 2 < indices.size(); i += 3) {
            tris.push_back({indices[i], indices[i + 1], indices[i + 2]});
        }

        write_triangulation(polygon, tris, args.output);
        std::cout << "seidel_polytri,vertices=" << polygon.size()
                  << ",triangles=" << tris.size()
                  << ",time_ms=" << elapsed_ms << "\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

