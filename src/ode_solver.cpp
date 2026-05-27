#include "ode_solver.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace af2d {

namespace {

constexpr double EPS = 1e-12;

// For Linear F(x) = slope*x + intercept, solve M'(t) = slope*M + intercept
// with M(t0) = x0. Return t1 such that M(t1) = x1.
double solve_linear_piece(double slope, double intercept, double t0, double x0, double x1) {
    if (std::abs(slope) < EPS) {
        // M(t) = x0 + intercept*(t - t0)
        if (std::abs(intercept) < EPS) {
            // F ~= 0, can't make progress
            return std::numeric_limits<double>::infinity();
        }
        return t0 + (x1 - x0) / intercept;
    }
    // M(t) = (x0 + intercept/slope) * exp(slope*(t-t0)) - intercept/slope
    // x1 + intercept/slope = (x0 + intercept/slope) * exp(slope*(t1-t0))
    double ratio = (x1 + intercept / slope) / (x0 + intercept / slope);
    if (ratio <= 0) {
        return std::numeric_limits<double>::infinity();
    }
    return t0 + std::log(ratio) / slope;
}

double eval_linear_solution(double slope, double intercept, double t0, double x0, double t) {
    if (std::abs(slope) < EPS) {
        return x0 + intercept * (t - t0);
    }
    return (x0 + intercept / slope) * std::exp(slope * (t - t0)) - intercept / slope;
}

// For SqrtQuadratic F(x) = sqrt((x-c)^2 + d^2), solve M'(t) = F(M(t))
// with M(t0) = x0. Uses M(t) = c + d*sinh(t - t0 + arcsinh((x0-c)/d))
double solve_sqrt_quad_piece(double c, double d, double t0, double x0, double x1) {
    if (std::abs(d) < EPS) {
        // F(x) = |x - c|. If x > c: M'(t) = M - c => M(t) = c + (x0-c)*exp(t-t0)
        // If x0 >= c (and x1 >= c since M is increasing):
        if (x0 >= c - EPS) {
            if (std::abs(x0 - c) < EPS) {
                // Starting exactly at c with F=0, can't make progress
                return std::numeric_limits<double>::infinity();
            }
            return t0 + std::log((x1 - c) / (x0 - c));
        }
        // x0 < c: M'(t) = c - M => M(t) = c - (c - x0)*exp(-(t-t0))
        // This approaches c but never reaches it, let alone x1 > c
        // In practice this shouldn't happen if LFS > 0
        return std::numeric_limits<double>::infinity();
    }
    // M(t) = c + d * sinh(t - t0 + arcsinh((x0 - c) / d))
    double asinh_arg = (x0 - c) / d;
    double offset = std::asinh(asinh_arg);
    // M(t1) = x1 => sinh(t1 - t0 + offset) = (x1 - c) / d
    double target_sinh = (x1 - c) / d;
    double t1 = t0 + std::asinh(target_sinh) - offset;
    return t1;
}

double eval_sqrt_quad_solution(double c, double d, double t0, double x0, double t) {
    if (std::abs(d) < EPS) {
        if (x0 >= c - EPS) {
            if (std::abs(x0 - c) < EPS) return c;
            return c + (x0 - c) * std::exp(t - t0);
        }
        return c - (c - x0) * std::exp(-(t - t0));
    }
    double offset = std::asinh((x0 - c) / d);
    return c + d * std::sinh(t - t0 + offset);
}

} // anonymous namespace

double ODESolution::evaluate(double t) const {
    if (pieces.empty()) return 0.0;

    for (const auto& piece : pieces) {
        if (t >= piece.t_start - EPS && t <= piece.t_end + EPS) {
            double tc = std::clamp(t, piece.t_start, piece.t_end);
            return std::visit([&](const auto& f) -> double {
                using T = std::decay_t<decltype(f)>;
                if constexpr (std::is_same_v<T, Linear>) {
                    return eval_linear_solution(f.slope, f.intercept, piece.t_start, piece.x_start, tc);
                } else {
                    return eval_sqrt_quad_solution(f.c, f.d, piece.t_start, piece.x_start, tc);
                }
            }, piece.func);
        }
    }

    // Beyond last piece
    return pieces.back().x_end;
}

ODESolution solve_ode(const LowerEnvelope& envelope, double length) {
    ODESolution result;

    if (envelope.pieces.empty() || length < EPS) {
        result.total_T = 0.0;
        return result;
    }

    double current_t = 0.0;
    double current_x = 0.0;

    for (const auto& env_piece : envelope.pieces) {
        // How far along x does this envelope piece cover?
        double piece_x_start = env_piece.x_start;
        double piece_x_end = env_piece.x_end;

        // We may not start at the beginning of this envelope piece
        if (current_x >= piece_x_end - EPS) continue;
        double x0 = std::max(current_x, piece_x_start);
        double x1 = std::min(piece_x_end, length);

        if (x1 <= x0 + EPS) continue;

        double t_end = std::visit([&](const auto& f) -> double {
            using T = std::decay_t<decltype(f)>;
            if constexpr (std::is_same_v<T, Linear>) {
                return solve_linear_piece(f.slope, f.intercept, current_t, x0, x1);
            } else {
                return solve_sqrt_quad_piece(f.c, f.d, current_t, x0, x1);
            }
        }, env_piece.func);

        if (std::isinf(t_end) || std::isnan(t_end)) {
            // Can't make progress; this shouldn't happen with valid LFS > 0
            break;
        }

        result.pieces.push_back({current_t, t_end, x0, x1, env_piece.func});
        current_t = t_end;
        current_x = x1;

        if (current_x >= length - EPS) break;
    }

    result.total_T = current_t;
    return result;
}

} // namespace af2d
