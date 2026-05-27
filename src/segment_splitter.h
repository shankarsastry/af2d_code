#pragma once

#include "pslg.h"
#include <string>
#include <vector>

namespace af2d {

struct SplitResult {
    std::vector<Point2D> all_vertices;   // original + split points
    std::vector<Segment> all_segments;   // sub-segments replacing originals
    std::vector<Point2D> holes;          // hole points passed through from PSLG
};

struct SplitParams {
    // Number of subsegments for the shortest reference segment.
    // Other segments get proportionally more: n_i = floor(n_star * T_i / T_min).
    //
    // Three ways to use:
    //   1. Leave at default (1): minimal splitting, subsegments sized by LFS ratio.
    //   2. Set directly:  params.n_star = 4;
    //   3. Compute from a minimum angle:
    //        auto result = compute_n_star_from_angle(20.0);
    //        params.n_star = result.n_star;
    //        // check result.message if needed
    int n_star = 1;
};

// Result from compute_n_star_from_angle.
struct NStarResult {
    int n_star;
    std::string message;  // non-empty if angle >= 30 or other issues
};

// Compute n* from the CDT non-encroachment condition (Lemma 4.12, Sastry 2021):
//   B*/A* + alpha/(A*(alpha-1)) + 2/A* <= 2*cos(theta*)
// where alpha = 1/(2*sin(theta*)), A* = n*/t_min_star - 1,
//       B* = n*/t_min_star + 1, t_min_star = 2*ln(2).
//
// For angles >= 30 degrees, the formula has a singularity (alpha -> 1)
// and cannot compute a meaningful n*. Returns n_star=1 with a message.
//
// Throws std::invalid_argument for angles <= 0 or >= 90.
NStarResult compute_n_star_from_angle(double theta_degrees);

// Split PSLG segments so subsegment lengths are proportional to local feature size.
SplitResult split_segments(const PSLG& pslg, const SplitParams& params = {});

} // namespace af2d
