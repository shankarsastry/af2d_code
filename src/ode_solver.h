#pragma once

#include "lower_envelope.h"
#include <vector>

namespace af2d {

struct PieceSolution {
    double t_start;
    double t_end;
    double x_start;
    double x_end;
    std::variant<SqrtQuadratic, Linear> func;  // the RHS F that was solved
};

struct ODESolution {
    std::vector<PieceSolution> pieces;
    double total_T = 0.0;
    double evaluate(double t) const;
};

// Solve M'(t) = F(M(t)) where F is given by the lower envelope.
// M(0) = 0, integrate until M reaches `length`.
ODESolution solve_ode(const LowerEnvelope& envelope, double length);

} // namespace af2d
