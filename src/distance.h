#pragma once

#include "pslg.h"
#include <variant>
#include <vector>

namespace af2d {

struct SqrtQuadratic {
    double c;  // center
    double d;  // perpendicular distance (>= 0)
    double evaluate(double x) const;
};

struct Linear {
    double slope;
    double intercept;
    double evaluate(double x) const;
};

struct DistancePiece {
    double x_start;
    double x_end;
    std::variant<SqrtQuadratic, Linear> func;
    double evaluate(double x) const;
};

struct DistanceFunction {
    std::vector<DistancePiece> pieces;
    double evaluate(double x) const;
};

// Distance from point at parameter x on segment seg_i to vertex vert_idx
DistanceFunction distance_to_vertex(const PSLG& pslg, int seg_i, int vert_idx);

// Distance from point at parameter x on segment seg_i to non-adjacent segment seg_j
DistanceFunction distance_to_segment(const PSLG& pslg, int seg_i, int seg_j);

// Distance to farthest endpoint of a segment: max(x, l-x) where l is segment length.
// This caps the local feature size per the paper's definition (Section 2.2, Fig. 5c).
DistanceFunction distance_to_farthest_endpoint(double length);

} // namespace af2d
