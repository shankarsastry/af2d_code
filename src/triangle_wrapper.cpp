#include "triangle_wrapper.h"

// Triangle library header (compiled as C++ by libigl mirror, no extern "C")
#include "triangle.h"

#include <cstring>
#include <stdexcept>

namespace af2d {

namespace {

// RAII wrapper for triangulateio's Triangle-allocated arrays
struct TriangulateIOGuard {
    triangulateio io{};

    TriangulateIOGuard() {
        std::memset(&io, 0, sizeof(io));
    }

    ~TriangulateIOGuard() {
        if (io.pointlist) free(io.pointlist);
        if (io.pointattributelist) free(io.pointattributelist);
        if (io.pointmarkerlist) free(io.pointmarkerlist);
        if (io.trianglelist) free(io.trianglelist);
        if (io.triangleattributelist) free(io.triangleattributelist);
        if (io.neighborlist) free(io.neighborlist);
        if (io.segmentlist) free(io.segmentlist);
        if (io.segmentmarkerlist) free(io.segmentmarkerlist);
        if (io.edgelist) free(io.edgelist);
        if (io.edgemarkerlist) free(io.edgemarkerlist);
        if (io.holelist) free(io.holelist);
        if (io.regionlist) free(io.regionlist);
    }

    TriangulateIOGuard(const TriangulateIOGuard&) = delete;
    TriangulateIOGuard& operator=(const TriangulateIOGuard&) = delete;
};

} // anonymous namespace

TriangleOutput triangulate_pslg(const SplitResult& split, const std::string& switches) {
    if (split.all_vertices.size() < 3) {
        throw std::invalid_argument("Need at least 3 vertices for triangulation");
    }

    // Set up input
    triangulateio in;
    std::memset(&in, 0, sizeof(in));

    int np = static_cast<int>(split.all_vertices.size());
    int ns = static_cast<int>(split.all_segments.size());

    std::vector<double> pointlist(np * 2);
    for (int i = 0; i < np; ++i) {
        pointlist[2 * i] = split.all_vertices[i].x;
        pointlist[2 * i + 1] = split.all_vertices[i].y;
    }

    std::vector<int> segmentlist(ns * 2);
    for (int i = 0; i < ns; ++i) {
        segmentlist[2 * i] = split.all_segments[i].p;
        segmentlist[2 * i + 1] = split.all_segments[i].q;
    }

    in.numberofpoints = np;
    in.pointlist = pointlist.data();
    in.numberofsegments = ns;
    in.segmentlist = segmentlist.data();
    int nh = static_cast<int>(split.holes.size());
    std::vector<double> holelist(nh * 2);
    for (int i = 0; i < nh; ++i) {
        holelist[2 * i] = split.holes[i].x;
        holelist[2 * i + 1] = split.holes[i].y;
    }
    in.numberofholes = nh;
    in.holelist = nh > 0 ? holelist.data() : nullptr;
    in.numberofregions = 0;

    // Set up output
    TriangulateIOGuard out;

    // Call Triangle
    // Need mutable string for Triangle's API
    std::string sw = switches;
    ::triangulate(sw.data(), &in, &out.io, nullptr);

    // Triangle copies in->holelist pointer to out->holelist.
    // Null it so TriangulateIOGuard doesn't free our stack data.
    out.io.holelist = nullptr;

    // Convert output
    TriangleOutput result;
    result.num_points = out.io.numberofpoints;
    result.num_triangles = out.io.numberoftriangles;
    result.num_edges = out.io.numberofedges;

    result.points.resize(result.num_points);
    for (int i = 0; i < result.num_points; ++i) {
        result.points[i] = {out.io.pointlist[2 * i], out.io.pointlist[2 * i + 1]};
    }

    result.triangles.resize(result.num_triangles);
    for (int i = 0; i < result.num_triangles; ++i) {
        result.triangles[i] = {
            out.io.trianglelist[3 * i],
            out.io.trianglelist[3 * i + 1],
            out.io.trianglelist[3 * i + 2]
        };
    }

    if (out.io.edgelist) {
        result.edges.resize(result.num_edges);
        for (int i = 0; i < result.num_edges; ++i) {
            result.edges[i] = {out.io.edgelist[2 * i], out.io.edgelist[2 * i + 1]};
        }
    }

    return result;
}

TriangleOutput triangulate(const PSLG& pslg, const SplitParams& params, const std::string& switches) {
    SplitResult split = split_segments(pslg, params);
    return triangulate_pslg(split, switches);
}

} // namespace af2d
