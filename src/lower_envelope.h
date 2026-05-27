#pragma once

#include "distance.h"
#include <vector>

namespace af2d {

struct EnvelopePiece {
    double x_start;
    double x_end;
    std::variant<SqrtQuadratic, Linear> func;
    double evaluate(double x) const;
};

struct LowerEnvelope {
    std::vector<EnvelopePiece> pieces;
    double evaluate(double x) const;
};

// Compute the pointwise minimum of multiple distance functions over [0, length]
LowerEnvelope compute_lower_envelope(const std::vector<DistanceFunction>& functions, double length);

} // namespace af2d
