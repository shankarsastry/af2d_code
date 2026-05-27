#include "distance.h"
#include <algorithm>
#include <cmath>

namespace af2d {

double SqrtQuadratic::evaluate(double x) const {
    return std::sqrt((x - c) * (x - c) + d * d);
}

double Linear::evaluate(double x) const {
    return slope * x + intercept;
}

double DistancePiece::evaluate(double x) const {
    return std::visit([x](const auto& f) { return f.evaluate(x); }, func);
}

double DistanceFunction::evaluate(double x) const {
    for (const auto& piece : pieces) {
        if (x >= piece.x_start - 1e-12 && x <= piece.x_end + 1e-12) {
            return piece.evaluate(x);
        }
    }
    // Fallback: evaluate nearest piece
    if (x <= pieces.front().x_start) return pieces.front().evaluate(pieces.front().x_start);
    return pieces.back().evaluate(pieces.back().x_end);
}

DistanceFunction distance_to_vertex(const PSLG& pslg, int seg_i, int vert_idx) {
    const auto& seg = pslg.segments[seg_i];
    const auto& pa = pslg.vertices[seg.p];
    const auto& pb = pslg.vertices[seg.q];
    const auto& v = pslg.vertices[vert_idx];

    Point2D dir = pb - pa;
    double l = dir.norm();
    Point2D t_hat = dir / l;
    Point2D n_hat = {-t_hat.y, t_hat.x};

    Point2D vp = v - pa;
    double a = vp.dot(t_hat);
    double b = std::abs(vp.dot(n_hat));

    DistanceFunction result;
    result.pieces.push_back({0.0, l, SqrtQuadratic{a, b}});
    return result;
}

DistanceFunction distance_to_segment(const PSLG& pslg, int seg_i, int seg_j) {
    const auto& si = pslg.segments[seg_i];
    const auto& pa = pslg.vertices[si.p];
    const auto& pb = pslg.vertices[si.q];

    const auto& sj = pslg.segments[seg_j];
    const auto& qa = pslg.vertices[sj.p];
    const auto& qb = pslg.vertices[sj.q];

    Point2D dir_i = pb - pa;
    double l_i = dir_i.norm();
    Point2D t_hat = dir_i / l_i;
    Point2D n_hat = {-t_hat.y, t_hat.x};

    // Project S_j endpoints into S_i's local frame
    Point2D qa_rel = qa - pa;
    Point2D qb_rel = qb - pa;
    double qa_t = qa_rel.dot(t_hat);
    double qa_n = qa_rel.dot(n_hat);
    double qb_t = qb_rel.dot(t_hat);
    double qb_n = qb_rel.dot(n_hat);

    // Direction of S_j in S_i's local frame
    double dt = qb_t - qa_t;
    double dn = qb_n - qa_n;
    double l_j = std::sqrt(dt * dt + dn * dn);

    // For each x on S_i, the closest point on S_j is found by projecting.
    // Parameterize S_j as Q(s) = qa + s*(qb-qa), s in [0,1].
    // In local frame: Q_t(s) = qa_t + s*dt, Q_n(s) = qa_n + s*dn
    // Distance^2 from (x,0) to Q(s): (x - qa_t - s*dt)^2 + (qa_n + s*dn)^2
    // Minimize over s: d/ds = 0 => s*(dt^2+dn^2) = (x-qa_t)*dt - qa_n*dn
    // s_opt(x) = ((x-qa_t)*dt - qa_n*dn) / (dt^2+dn^2)

    // The distance function is piecewise:
    // - When s_opt < 0: distance to near endpoint qa
    // - When s_opt in [0,1]: perpendicular distance (which is linear or sqrt-quad)
    // - When s_opt > 1: distance to far endpoint qb

    double denom = dt * dt + dn * dn; // = l_j^2

    DistanceFunction result;

    if (denom < 1e-24) {
        // Degenerate segment j (zero length): treat as vertex
        result.pieces.push_back({0.0, l_i, SqrtQuadratic{qa_t, std::abs(qa_n)}});
        return result;
    }

    // Find x values where s_opt = 0 and s_opt = 1
    // s_opt(x) = ((x-qa_t)*dt - qa_n*dn) / denom
    // s_opt = 0 => x = qa_t + qa_n*dn/dt ... but better:
    // s_opt(x) = 0 => (x-qa_t)*dt = qa_n*dn => x = qa_t + qa_n*dn/dt
    // But this only makes sense when dt != 0. Let's use the general form:
    // s_opt(x) * denom = (x - qa_t)*dt - qa_n*dn
    // s_opt = 0 at x0 = qa_t + qa_n*dn/dt  (if dt != 0)
    // s_opt = 1 at x1 = qa_t + (denom + qa_n*dn)/dt = qa_t + dt + qa_n*dn/dt (if dt != 0)

    // More carefully: s_opt(x) = ((x - qa_t)*dt - qa_n*dn) / denom
    // s_opt(x) = 0 => x*dt = qa_t*dt + qa_n*dn => x = qa_t + qa_n*dn/dt
    // s_opt(x) = 1 => (x-qa_t)*dt - qa_n*dn = denom => x = qa_t + (denom + qa_n*dn)/dt

    // If dt == 0 (segments are perpendicular in tangent direction), s_opt is constant.

    double x_s0, x_s1; // x where s_opt = 0 and s_opt = 1

    if (std::abs(dt) < 1e-15) {
        // s_opt doesn't depend on x: s_opt = -qa_n*dn / denom
        double s_const = -qa_n * dn / denom;
        if (s_const <= 0.0) {
            // Always closest to qa
            result.pieces.push_back({0.0, l_i, SqrtQuadratic{qa_t, std::abs(qa_n)}});
            return result;
        } else if (s_const >= 1.0) {
            // Always closest to qb
            result.pieces.push_back({0.0, l_i, SqrtQuadratic{qb_t, std::abs(qb_n)}});
            return result;
        } else {
            // Closest to interior point. Perpendicular distance is constant.
            // Q(s_const) in local frame: (qa_t + s_const*dt, qa_n + s_const*dn) = (qa_t, qa_n + s_const*dn)
            // since dt=0. Distance from (x,0) to (qa_t, qa_n+s_const*dn) = sqrt((x-qa_t)^2 + (qa_n+s_const*dn)^2)
            // But wait, if dt=0 and the closest point is interior, the perpendicular from S_i to S_j
            // is at a fixed n-coordinate. The distance is |n_perp|... Let me reconsider.
            // Actually with dt=0, segments are "parallel" in a sense (S_j is purely normal-direction).
            // The closest point on S_j to (x,0) has s_opt = -qa_n*dn/denom (constant).
            // The point is at (qa_t, qa_n + s_const*dn).
            double n_closest = qa_n + s_const * dn;
            result.pieces.push_back({0.0, l_i, SqrtQuadratic{qa_t, std::abs(n_closest)}});
            return result;
        }
    }

    x_s0 = qa_t + qa_n * dn / dt;
    x_s1 = qa_t + (denom + qa_n * dn) / dt;

    // Ensure x_s0 < x_s1 (swap if dt < 0)
    if (x_s0 > x_s1) std::swap(x_s0, x_s1);

    // Build piecewise function over [0, l_i]
    // Region 1: x < x_s0 => closest to one endpoint
    // Region 2: x_s0 <= x <= x_s1 => closest to interior of S_j
    // Region 3: x > x_s1 => closest to other endpoint

    // Determine which endpoint is "near" (s=0 side) vs "far" (s=1 side)
    // When dt > 0, x < x_s0 => s_opt < 0 => closest to qa, x > x_s1 => s_opt > 1 => closest to qb
    // When dt < 0, it's swapped, but we already swapped x_s0/x_s1
    int near_vert, far_vert;
    double near_t, near_n, far_t, far_n;
    if (dt > 0) {
        near_t = qa_t; near_n = qa_n;
        far_t = qb_t; far_n = qb_n;
    } else {
        near_t = qb_t; near_n = qb_n;
        far_t = qa_t; far_n = qa_n;
    }

    // In the interior region, the closest point on S_j traces along S_j.
    // The distance is: let's compute it.
    // The perpendicular distance from the line containing S_i to the line containing S_j:
    // Actually the distance from (x, 0) to Q(s_opt(x)) where s_opt is in [0,1].
    // Let me compute this differently for the interior region.
    //
    // Q(s) = (qa_t + s*dt, qa_n + s*dn), and s_opt(x) = ((x-qa_t)*dt - qa_n*dn)/denom
    // Distance^2 = (x - qa_t - s*dt)^2 + (qa_n + s*dn)^2 at s=s_opt
    // After substitution and simplification:
    // D^2 = ((x-qa_t)*dn - qa_n*(-dt))^2 / denom ... let me just compute it properly.
    //
    // With s = s_opt: x - qa_t - s*dt = x - qa_t - ((x-qa_t)*dt - qa_n*dn)*dt/denom
    //   = (x-qa_t)*(1 - dt^2/denom) + qa_n*dn*dt/denom
    //   = (x-qa_t)*dn^2/denom + qa_n*dn*dt/denom
    //   = dn*((x-qa_t)*dn + qa_n*dt)/denom
    //
    // qa_n + s*dn = qa_n + ((x-qa_t)*dt - qa_n*dn)*dn/denom
    //   = (qa_n*denom + (x-qa_t)*dt*dn - qa_n*dn^2)/denom
    //   = (qa_n*dt^2 + (x-qa_t)*dt*dn)/denom
    //   = dt*(qa_n*dt + (x-qa_t)*dn)/denom
    //
    // D^2 = (dn^2 + dt^2)*((x-qa_t)*dn + qa_n*dt)^2 / denom^2
    //      = denom * ((x-qa_t)*dn + qa_n*dt)^2 / denom^2
    //      = ((x-qa_t)*dn + qa_n*dt)^2 / denom
    //
    // D = |((x-qa_t)*dn + qa_n*dt)| / sqrt(denom)
    //   = |dn*(x-qa_t) + qa_n*dt| / l_j
    //
    // This is a linear function of x! slope = |dn|/l_j, but with absolute value.
    // D(x) = |dn*x - dn*qa_t + qa_n*dt| / l_j
    //       = |dn*x + (qa_n*dt - dn*qa_t)| / l_j

    double perp_slope = dn / l_j;  // can be negative
    double perp_const = (qa_n * dt - dn * qa_t) / l_j;

    // D(x) = |perp_slope * x + perp_const|
    // For a non-degenerate non-adjacent segment, this should be positive.
    // If segments don't intersect, sign is constant => D(x) = ±(perp_slope*x + perp_const)

    // Check sign at midpoint of the interior region
    double x_mid_interior = std::clamp((x_s0 + x_s1) / 2.0, 0.0, l_i);
    double val_at_mid = perp_slope * x_mid_interior + perp_const;
    double sign = (val_at_mid >= 0) ? 1.0 : -1.0;
    double eff_slope = sign * perp_slope;
    double eff_intercept = sign * perp_const;

    // Now build the pieces
    double domain_start = 0.0;
    double domain_end = l_i;

    // Clamp boundary points to domain
    double b0 = std::clamp(x_s0, domain_start, domain_end);
    double b1 = std::clamp(x_s1, domain_start, domain_end);

    // Region 1: [0, b0] -> sqrt-quadratic to near endpoint
    if (b0 > domain_start + 1e-15) {
        result.pieces.push_back({domain_start, b0, SqrtQuadratic{near_t, std::abs(near_n)}});
    }

    // Region 2: [b0, b1] -> linear perpendicular distance
    if (b1 > b0 + 1e-15) {
        result.pieces.push_back({b0, b1, Linear{eff_slope, eff_intercept}});
    }

    // Region 3: [b1, l_i] -> sqrt-quadratic to far endpoint
    if (domain_end > b1 + 1e-15) {
        result.pieces.push_back({b1, domain_end, SqrtQuadratic{far_t, std::abs(far_n)}});
    }

    // If no pieces were created (degenerate case), add a full-domain distance to nearest endpoint
    if (result.pieces.empty()) {
        // Both boundary points outside domain - determine which endpoint is relevant
        if (x_s1 <= domain_start) {
            // All of S_i maps to s>1 region (far endpoint)
            result.pieces.push_back({domain_start, domain_end, SqrtQuadratic{far_t, std::abs(far_n)}});
        } else if (x_s0 >= domain_end) {
            // All of S_i maps to s<0 region (near endpoint)
            result.pieces.push_back({domain_start, domain_end, SqrtQuadratic{near_t, std::abs(near_n)}});
        } else {
            // Shouldn't happen, but fallback
            result.pieces.push_back({domain_start, domain_end, Linear{eff_slope, eff_intercept}});
        }
    }

    return result;
}

DistanceFunction distance_to_farthest_endpoint(double length) {
    double half = length / 2.0;
    DistanceFunction result;
    // [0, l/2]: farthest endpoint is q at distance l-x => slope=-1, intercept=l
    result.pieces.push_back({0.0, half, Linear{-1.0, length}});
    // [l/2, l]: farthest endpoint is p at distance x => slope=1, intercept=0
    result.pieces.push_back({half, length, Linear{1.0, 0.0}});
    return result;
}

} // namespace af2d
