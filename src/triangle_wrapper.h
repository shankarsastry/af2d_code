#pragma once

#include "segment_splitter.h"
#include <array>
#include <string>
#include <vector>

namespace af2d {

struct TriangleOutput {
    std::vector<Point2D> points;
    std::vector<std::array<int, 3>> triangles;
    std::vector<Segment> edges;
    int num_points = 0;
    int num_triangles = 0;
    int num_edges = 0;
};

// Triangulate a split PSLG using Shewchuk's Triangle library.
// switches: Triangle command-line switches (default "pzQ" = PSLG, zero-indexed, quiet)
TriangleOutput triangulate_pslg(const SplitResult& split,
                                 const std::string& switches = "pzQ");

// Convenience: triangulate directly from a PSLG (splits first)
TriangleOutput triangulate(const PSLG& pslg, const SplitParams& params = {},
                           const std::string& switches = "pzQ");

} // namespace af2d
