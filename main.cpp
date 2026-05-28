#include "io.h"
#include "pslg.h"
#include "segment_splitter.h"
#include "triangle_wrapper.h"
#include <cstring>
#include <iostream>
#include <string>

using namespace af2d;

static void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <input.poly|input.node> [-q<angle>] [-n<n_star>] [-j<threads>]\n"
              << "\n"
              << "  <input>       Input file (.poly or .node)\n"
              << "  -q<angle>     Minimum angle in degrees (computes n* if < 30)\n"
              << "  -n<n_star>    Explicit n* value (overrides -q)\n"
              << "  -t            Enable quadtree acceleration\n"
              << "  -j<threads>   Number of OpenMP threads (0 = default)\n"
              << "\n"
              << "Output: <basename>.1.node, <basename>.1.ele, <basename>.1.edge\n";
}

// Strip the extension (.poly or .node) from a path to get the basename.
static std::string strip_extension(const std::string& path) {
    auto dot = path.rfind('.');
    if (dot == std::string::npos) return path;
    return path.substr(0, dot);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string input_path;
    double q_angle = -1.0;
    int explicit_n_star = -1;
    bool use_quadtree = false;
    int num_threads = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-t") {
            use_quadtree = true;
        } else if (arg.rfind("-q", 0) == 0) {
            q_angle = std::stod(arg.substr(2));
        } else if (arg.rfind("-n", 0) == 0) {
            explicit_n_star = std::stoi(arg.substr(2));
        } else if (arg.rfind("-j", 0) == 0) {
            num_threads = std::stoi(arg.substr(2));
        } else if (arg[0] == '-') {
            std::cerr << "Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        } else {
            input_path = arg;
        }
    }

    if (input_path.empty()) {
        std::cerr << "Error: no input file specified.\n";
        print_usage(argv[0]);
        return 1;
    }

    // Read input
    PSLG pslg;
    bool is_poly = (input_path.size() >= 5 &&
                    input_path.substr(input_path.size() - 5) == ".poly");
    if (is_poly) {
        pslg = read_poly(input_path);
    } else {
        pslg = read_node(input_path);
    }

    pslg.validate();

    // Determine n*
    SplitParams params;
    if (explicit_n_star > 0) {
        params.n_star = explicit_n_star;
    } else if (q_angle > 0.0) {
        auto nsr = compute_n_star_from_angle(q_angle);
        params.n_star = nsr.n_star;
        if (!nsr.message.empty()) {
            std::cerr << "Warning: " << nsr.message << "\n";
        }
    }
    // else default n_star = 1

    params.use_quadtree = use_quadtree;
    params.num_threads = num_threads;

    // Split and triangulate
    auto split = split_segments(pslg, params);

    // Build Triangle switches: always use pzQe (PSLG, zero-indexed, quiet, edges).
    // If -q was given, pass it to Triangle for minimum angle enforcement.
    std::string switches = "pzQe";
    if (q_angle > 0.0) {
        switches += "q" + std::to_string(q_angle);
    }
    auto out = triangulate_pslg(split, switches);

    // Write output
    std::string basename = strip_extension(input_path) + ".1";
    write_node(basename + ".node", out);
    write_ele(basename + ".ele", out);
    write_edge(basename + ".edge", out);

    // Summary
    std::cout << "Input:     " << pslg.vertices.size() << " vertices, "
              << pslg.segments.size() << " segments, "
              << pslg.holes.size() << " holes\n";
    std::cout << "n*:        " << params.n_star << "\n";
    std::cout << "Output:    " << out.num_points << " vertices, "
              << out.num_triangles << " triangles, "
              << out.num_edges << " edges\n";
    std::cout << "Files:     " << basename << ".node\n"
              << "           " << basename << ".ele\n"
              << "           " << basename << ".edge\n";

    return 0;
}
