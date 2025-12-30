#include "polygon_io.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <chrono>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using Clock = std::chrono::high_resolution_clock;

namespace baselines {
namespace cgal_runner {

struct FaceInfo2 {
    int nesting_level;
    bool in_domain() const { return nesting_level % 2 == 1; }
};

} // namespace cgal_runner
} // namespace baselines

namespace {

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, K>;
using Fb_base = CGAL::Constrained_triangulation_face_base_2<K>;
using Fb = CGAL::Triangulation_face_base_with_info_2<
    baselines::cgal_runner::FaceInfo2,
    K,
    Fb_base>;
using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using Itag = CGAL::Exact_predicates_tag;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>;
using CGALPoint = K::Point_2;
using Edge = CDT::Edge;
using Face_handle = CDT::Face_handle;

void mark_domains(CDT& cdt,
                  Face_handle start,
                  int index,
                  std::list<Edge>& border) {
    if (start->info().nesting_level != -1) {
        return;
    }
    std::list<Face_handle> queue;
    queue.push_back(start);
    while (!queue.empty()) {
        Face_handle fh = queue.front();
        queue.pop_front();
        if (fh->info().nesting_level == -1) {
            fh->info().nesting_level = index;
            for (int i = 0; i < 3; i++) {
                Edge e(fh, i);
                Face_handle n = fh->neighbor(i);
                if (n->info().nesting_level == -1) {
                    if (cdt.is_constrained(e)) {
                        border.push_back(e);
                    } else {
                        queue.push_back(n);
                    }
                }
            }
        }
    }
}

void mark_domains(CDT& cdt) {
    for (auto fit = cdt.all_faces_begin(); fit != cdt.all_faces_end(); ++fit) {
        fit->info().nesting_level = -1;
    }

    std::list<Edge> border;
    mark_domains(cdt, cdt.infinite_face(), 0, border);
    while (!border.empty()) {
        Edge e = border.front();
        border.pop_front();
        Face_handle n = e.first->neighbor(e.second);
        if (n->info().nesting_level == -1) {
            mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
        }
    }
}

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

        CDT cdt;
        const std::size_t invalid = std::numeric_limits<std::size_t>::max();

        std::vector<CGALPoint> points;
        points.reserve(polygon.size());
        for (const auto& p : polygon) {
            points.emplace_back(p.x, p.y);
        }

        std::vector<CDT::Vertex_handle> handles;
        handles.reserve(points.size());
        for (std::size_t i = 0; i < points.size(); ++i) {
            auto vh = cdt.insert(points[i]);
            vh->info() = i;
            handles.push_back(vh);
        }

        for (std::size_t i = 0; i < handles.size(); ++i) {
            cdt.insert_constraint(handles[i], handles[(i + 1) % handles.size()]);
        }

        const auto start = Clock::now();
        mark_domains(cdt);
        const auto stop = Clock::now();
        const double elapsed_ms =
            std::chrono::duration<double, std::milli>(stop - start).count();

        // Assign indices to any new vertices (should be none for simple polygons)
        Polygon output_vertices = polygon;
        for (auto vit = cdt.finite_vertices_begin();
             vit != cdt.finite_vertices_end(); ++vit) {
            if (vit->info() == invalid) {
                const auto idx = output_vertices.size();
                vit->info() = idx;
                output_vertices.push_back({vit->point().x(), vit->point().y()});
            }
        }

        Triangles tris;
        for (auto fit = cdt.finite_faces_begin();
             fit != cdt.finite_faces_end(); ++fit) {
            if (!fit->info().in_domain()) continue;
            tris.push_back({
                static_cast<uint32_t>(fit->vertex(0)->info()),
                static_cast<uint32_t>(fit->vertex(1)->info()),
                static_cast<uint32_t>(fit->vertex(2)->info())
            });
        }

        write_triangulation(output_vertices, tris, args.output);
        std::cout << "cgal_cdt,vertices=" << polygon.size()
                  << ",triangles=" << tris.size()
                  << ",time_ms=" << elapsed_ms << "\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

