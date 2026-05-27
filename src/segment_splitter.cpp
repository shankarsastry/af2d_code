#include "segment_splitter.h"
#include "distance.h"
#include "lower_envelope.h"
#include "ode_solver.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace af2d {

namespace {

constexpr double T_MIN_STAR = 2.0 * 0.693147180559945309;  // 2*ln(2) ≈ 1.3863
constexpr double PI = 3.14159265358979323846;

struct SegmentODE {
    ODESolution ode;
    double length;
    bool has_features;  // whether non-adjacent features exist
};

SegmentODE compute_segment_ode(const PSLG& pslg, int seg_i) {
    SegmentODE result;
    result.length = pslg.segment_length(seg_i);
    result.has_features = false;

    std::vector<DistanceFunction> dist_funcs;

    for (int v = 0; v < static_cast<int>(pslg.vertices.size()); ++v) {
        if (!pslg.is_adjacent_vertex(seg_i, v)) {
            dist_funcs.push_back(distance_to_vertex(pslg, seg_i, v));
        }
    }

    for (int seg_j = 0; seg_j < static_cast<int>(pslg.segments.size()); ++seg_j) {
        if (seg_j != seg_i && !pslg.is_adjacent(seg_i, seg_j)) {
            dist_funcs.push_back(distance_to_segment(pslg, seg_i, seg_j));
        }
    }

    // Always include distance to farthest own endpoint: max(x, l-x).
    // This caps lfs per the paper (Section 2.2, Fig. 5c).
    dist_funcs.push_back(distance_to_farthest_endpoint(result.length));

    result.has_features = true;

    LowerEnvelope envelope = compute_lower_envelope(dist_funcs, result.length);
    if (envelope.pieces.empty()) {
        result.has_features = false;
        return result;
    }

    result.ode = solve_ode(envelope, result.length);
    if (result.ode.total_T < 1e-15 || result.ode.pieces.empty()) {
        result.has_features = false;
    }

    return result;
}

// Evaluate the CDT non-encroachment LHS for a given n_star.
// Lemma 4.12: B*/A* + alpha/(A*(alpha-1)) + 2/A* <= 2*cos(theta*)
double encroachment_lhs(int n_star, double theta_rad) {
    double alpha = 1.0 / (2.0 * std::sin(theta_rad));
    double a_star = static_cast<double>(n_star) / T_MIN_STAR - 1.0;
    double b_star = static_cast<double>(n_star) / T_MIN_STAR + 1.0;

    if (a_star <= 0.0) {
        return std::numeric_limits<double>::max();
    }

    double lhs = b_star / a_star + 2.0 / a_star;

    double denom = alpha - 1.0;
    if (std::abs(denom) > 1e-10) {
        lhs += alpha / (a_star * denom);
    } else {
        return std::numeric_limits<double>::max();
    }

    return lhs;
}

} // anonymous namespace

NStarResult compute_n_star_from_angle(double theta_degrees) {
    if (theta_degrees <= 0.0 || theta_degrees >= 90.0) {
        throw std::invalid_argument(
            "min_angle_degrees must be in (0, 90)");
    }

    if (theta_degrees >= 30.0) {
        return {1, "min_angle >= 30 degrees: the non-encroachment formula has a "
                   "singularity at alpha=1 (theta=30). n_star cannot be computed "
                   "from the angle; defaulting to 1."};
    }

    double theta_rad = theta_degrees * PI / 180.0;
    double rhs = 2.0 * std::cos(theta_rad);

    for (int n = 2; n <= 10000; ++n) {
        if (encroachment_lhs(n, theta_rad) <= rhs) {
            return {n, ""};
        }
    }

    return {1, "Could not find a feasible n_star within 10000 iterations; "
               "defaulting to 1."};
}

SplitResult split_segments(const PSLG& pslg, const SplitParams& params) {
    int num_segs = static_cast<int>(pslg.segments.size());
    int n_star = std::max(1, params.n_star);

    // Pass 1: compute ODE solutions for all segments
    std::vector<SegmentODE> seg_odes(num_segs);
    #pragma omp parallel for schedule(dynamic) if(num_segs > 1)
    for (int i = 0; i < num_segs; ++i) {
        seg_odes[i] = compute_segment_ode(pslg, i);
    }

    // Find the minimum total_T among segments that have non-adjacent features
    double min_T = std::numeric_limits<double>::max();
    for (int i = 0; i < num_segs; ++i) {
        if (seg_odes[i].has_features && seg_odes[i].ode.total_T > 1e-15) {
            min_T = std::min(min_T, seg_odes[i].ode.total_T);
        }
    }

    // Pass 2: split each segment
    // For segment i with total_T = T_i:
    //   n_i = floor(n_star * T_i / min_T)
    // This gives n_i = n_star for the shortest segment, proportionally more for longer ones.
    SplitResult result;
    result.all_vertices = pslg.vertices;
    result.holes = pslg.holes;
    int next_vertex_idx = static_cast<int>(result.all_vertices.size());

    for (int seg_i = 0; seg_i < num_segs; ++seg_i) {
        const auto& seg_ode = seg_odes[seg_i];

        if (!seg_ode.has_features) {
            result.all_segments.push_back(pslg.segments[seg_i]);
            continue;
        }

        double T_i = seg_ode.ode.total_T;
        int n_i;
        if (min_T < std::numeric_limits<double>::max() && min_T > 1e-15) {
            n_i = std::max(1, static_cast<int>(std::floor(
                static_cast<double>(n_star) * T_i / min_T)));
        } else {
            n_i = 1;
        }

        double h = T_i / n_i;

        std::vector<int> chain_vertices;
        chain_vertices.push_back(pslg.segments[seg_i].p);

        for (int k = 1; k < n_i; ++k) {
            double t_k = k * h;
            double x_k = seg_ode.ode.evaluate(t_k);
            x_k = std::clamp(x_k, 0.0, seg_ode.length);

            Point2D pt = pslg.point_on_segment(seg_i, x_k);
            result.all_vertices.push_back(pt);
            chain_vertices.push_back(next_vertex_idx++);
        }

        chain_vertices.push_back(pslg.segments[seg_i].q);

        for (size_t k = 0; k + 1 < chain_vertices.size(); ++k) {
            result.all_segments.push_back({chain_vertices[k], chain_vertices[k + 1]});
        }
    }

    return result;
}

} // namespace af2d
