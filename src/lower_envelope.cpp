#include "lower_envelope.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace af2d {

double EnvelopePiece::evaluate(double x) const {
    return std::visit([x](const auto& f) { return f.evaluate(x); }, func);
}

double LowerEnvelope::evaluate(double x) const {
    for (const auto& piece : pieces) {
        if (x >= piece.x_start - 1e-12 && x <= piece.x_end + 1e-12) {
            return piece.evaluate(x);
        }
    }
    if (pieces.empty()) return std::numeric_limits<double>::max();
    if (x <= pieces.front().x_start) return pieces.front().evaluate(pieces.front().x_start);
    return pieces.back().evaluate(pieces.back().x_end);
}

namespace {

constexpr double EPS = 1e-10;

// Find intersections between two function pieces on a given interval
std::vector<double> find_intersections(
    const std::variant<SqrtQuadratic, Linear>& f1,
    const std::variant<SqrtQuadratic, Linear>& f2,
    double x_lo, double x_hi)
{
    std::vector<double> roots;

    auto is_in_range = [&](double x) {
        return x > x_lo + EPS && x < x_hi - EPS;
    };

    if (auto* sq1 = std::get_if<SqrtQuadratic>(&f1)) {
        if (auto* sq2 = std::get_if<SqrtQuadratic>(&f2)) {
            // sqrt((x-c1)^2 + d1^2) = sqrt((x-c2)^2 + d2^2)
            // => (x-c1)^2 + d1^2 = (x-c2)^2 + d2^2
            // => x^2 - 2*c1*x + c1^2 + d1^2 = x^2 - 2*c2*x + c2^2 + d2^2
            // => 2*(c2-c1)*x = c2^2 - c1^2 + d2^2 - d1^2
            double dc = sq2->c - sq1->c;
            if (std::abs(dc) > EPS) {
                double x_cross = (sq2->c * sq2->c - sq1->c * sq1->c +
                                  sq2->d * sq2->d - sq1->d * sq1->d) / (2.0 * dc);
                if (is_in_range(x_cross)) roots.push_back(x_cross);
            }
            // If dc ~= 0 and d1 ~= d2, functions are identical (no crossing)
        } else {
            auto* lin = std::get_if<Linear>(&f2);
            // sqrt((x-c)^2 + d^2) = slope*x + intercept
            // Squaring: (x-c)^2 + d^2 = (slope*x + intercept)^2
            // x^2 - 2cx + c^2 + d^2 = slope^2*x^2 + 2*slope*intercept*x + intercept^2
            // (1-slope^2)*x^2 + (-2c - 2*slope*intercept)*x + (c^2+d^2-intercept^2) = 0
            double A = 1.0 - lin->slope * lin->slope;
            double B = -2.0 * sq1->c - 2.0 * lin->slope * lin->intercept;
            double C = sq1->c * sq1->c + sq1->d * sq1->d - lin->intercept * lin->intercept;
            double disc = B * B - 4.0 * A * C;
            if (disc >= 0) {
                double sq_disc = std::sqrt(std::max(0.0, disc));
                if (std::abs(A) > EPS) {
                    double x1 = (-B + sq_disc) / (2.0 * A);
                    double x2 = (-B - sq_disc) / (2.0 * A);
                    // Verify each root (squaring can introduce spurious roots)
                    for (double xr : {x1, x2}) {
                        if (is_in_range(xr)) {
                            double lhs = sq1->evaluate(xr);
                            double rhs = lin->evaluate(xr);
                            if (std::abs(lhs - rhs) < 1e-6 * std::max(1.0, std::abs(lhs))) {
                                roots.push_back(xr);
                            }
                        }
                    }
                } else {
                    // Linear in x: Bx + C = 0
                    if (std::abs(B) > EPS) {
                        double xr = -C / B;
                        if (is_in_range(xr)) {
                            double lhs = sq1->evaluate(xr);
                            double rhs = lin->evaluate(xr);
                            if (std::abs(lhs - rhs) < 1e-6 * std::max(1.0, std::abs(lhs))) {
                                roots.push_back(xr);
                            }
                        }
                    }
                }
            }
        }
    } else {
        auto* lin1 = std::get_if<Linear>(&f1);
        if (auto* sq2 = std::get_if<SqrtQuadratic>(&f2)) {
            // Reuse sqrt-quad vs linear
            return find_intersections(f2, f1, x_lo, x_hi);
        } else {
            auto* lin2 = std::get_if<Linear>(&f2);
            // slope1*x + int1 = slope2*x + int2
            double ds = lin1->slope - lin2->slope;
            if (std::abs(ds) > EPS) {
                double xr = (lin2->intercept - lin1->intercept) / ds;
                if (is_in_range(xr)) roots.push_back(xr);
            }
        }
    }

    return roots;
}

struct EvalFunc {
    const std::vector<DistanceFunction>* functions;
    double evaluate(int func_idx, double x) const {
        return (*functions)[func_idx].evaluate(x);
    }
};

} // anonymous namespace

LowerEnvelope compute_lower_envelope(const std::vector<DistanceFunction>& functions, double length) {
    LowerEnvelope result;
    if (functions.empty()) return result;

    if (functions.size() == 1) {
        for (const auto& piece : functions[0].pieces) {
            result.pieces.push_back({piece.x_start, piece.x_end, piece.func});
        }
        return result;
    }

    // Step 1: Collect all breakpoints
    std::set<double> breakpoints_set;
    breakpoints_set.insert(0.0);
    breakpoints_set.insert(length);

    for (const auto& df : functions) {
        for (const auto& piece : df.pieces) {
            breakpoints_set.insert(piece.x_start);
            breakpoints_set.insert(piece.x_end);
        }
    }

    // Step 2: Find pairwise intersections
    // For each pair of functions, find where they cross within the domain
    for (size_t i = 0; i < functions.size(); ++i) {
        for (size_t j = i + 1; j < functions.size(); ++j) {
            // Check each pair of pieces
            for (const auto& pi : functions[i].pieces) {
                for (const auto& pj : functions[j].pieces) {
                    double lo = std::max(pi.x_start, pj.x_start);
                    double hi = std::min(pi.x_end, pj.x_end);
                    if (hi <= lo + EPS) continue;

                    auto crossings = find_intersections(pi.func, pj.func, lo, hi);
                    for (double x : crossings) {
                        breakpoints_set.insert(x);
                    }
                }
            }
        }
    }

    // Step 3: Sort breakpoints
    std::vector<double> breakpoints(breakpoints_set.begin(), breakpoints_set.end());
    std::sort(breakpoints.begin(), breakpoints.end());

    // Remove breakpoints outside [0, length]
    breakpoints.erase(
        std::remove_if(breakpoints.begin(), breakpoints.end(),
                        [length](double x) { return x < -EPS || x > length + EPS; }),
        breakpoints.end());

    // Step 4: For each sub-interval, evaluate all functions at midpoint and pick minimum
    for (size_t k = 0; k + 1 < breakpoints.size(); ++k) {
        double x_lo = breakpoints[k];
        double x_hi = breakpoints[k + 1];
        if (x_hi - x_lo < EPS) continue;

        double x_mid = (x_lo + x_hi) / 2.0;

        double best_val = std::numeric_limits<double>::max();
        int best_idx = -1;

        for (size_t i = 0; i < functions.size(); ++i) {
            double val = functions[i].evaluate(x_mid);
            if (val < best_val - EPS) {
                best_val = val;
                best_idx = static_cast<int>(i);
            }
        }

        if (best_idx < 0) continue;

        // Find which piece of the winning function covers this interval
        const auto& winner = functions[best_idx];
        for (const auto& piece : winner.pieces) {
            if (x_mid >= piece.x_start - EPS && x_mid <= piece.x_end + EPS) {
                result.pieces.push_back({x_lo, x_hi, piece.func});
                break;
            }
        }
    }

    // Step 5: Merge adjacent intervals with same function type and parameters
    if (result.pieces.size() <= 1) return result;

    std::vector<EnvelopePiece> merged;
    merged.push_back(result.pieces[0]);

    for (size_t i = 1; i < result.pieces.size(); ++i) {
        auto& prev = merged.back();
        const auto& curr = result.pieces[i];

        bool can_merge = false;
        if (prev.func.index() == curr.func.index()) {
            if (auto* sq1 = std::get_if<SqrtQuadratic>(&prev.func)) {
                auto* sq2 = std::get_if<SqrtQuadratic>(&curr.func);
                can_merge = std::abs(sq1->c - sq2->c) < EPS && std::abs(sq1->d - sq2->d) < EPS;
            } else {
                auto* lin1 = std::get_if<Linear>(&prev.func);
                auto* lin2 = std::get_if<Linear>(&curr.func);
                can_merge = std::abs(lin1->slope - lin2->slope) < EPS &&
                            std::abs(lin1->intercept - lin2->intercept) < EPS;
            }
        }

        if (can_merge && std::abs(prev.x_end - curr.x_start) < EPS) {
            prev.x_end = curr.x_end;
        } else {
            merged.push_back(curr);
        }
    }

    result.pieces = std::move(merged);
    return result;
}

} // namespace af2d
