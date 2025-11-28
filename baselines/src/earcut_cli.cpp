#include "polygon_io.hpp"

#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "mapbox/earcut.hpp"

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

} // namespace

int main(int argc, char** argv) {
    using namespace baselines;
    try {
        const auto args = parse_args(argc, argv);
        const auto polygon = read_polygon(args.input);

        // Earcut input (single contour, no holes)
        std::vector<std::vector<std::array<double, 2>>> poly_data(1);
        poly_data[0].reserve(polygon.size());
        for (const auto& p : polygon) {
            poly_data[0].push_back({p.x, p.y});
        }

        const auto start = Clock::now();
        auto indices = mapbox::earcut<uint32_t>(poly_data);
        const auto stop = Clock::now();
        const double elapsed_ms =
            std::chrono::duration<double, std::milli>(stop - start).count();

        Triangles tris;
        tris.reserve(indices.size() / 3);
        for (size_t i = 0; i + 2 < indices.size(); i += 3) {
            tris.push_back({indices[i], indices[i + 1], indices[i + 2]});
        }

        write_triangulation(polygon, tris, args.output);
        std::cout << "earcut,vertices=" << polygon.size()
                  << ",triangles=" << tris.size()
                  << ",time_ms=" << elapsed_ms << "\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

