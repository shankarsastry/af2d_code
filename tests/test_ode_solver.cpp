#include <gtest/gtest.h>
#include "ode_solver.h"
#include <cmath>

using namespace af2d;

TEST(ODESolver, ConstantF) {
    // F(x) = c (constant), represented as Linear{0, c}
    // M'(t) = c => M(t) = c*t
    // total_T = length / c
    double c = 2.0;
    double length = 4.0;
    LowerEnvelope env;
    env.pieces.push_back({0.0, length, Linear{0.0, c}});

    auto sol = solve_ode(env, length);
    EXPECT_NEAR(sol.total_T, length / c, 1e-10);
    EXPECT_NEAR(sol.evaluate(0.0), 0.0, 1e-10);
    EXPECT_NEAR(sol.evaluate(sol.total_T), length, 1e-10);
}

TEST(ODESolver, LinearWithSlope) {
    // F(x) = x + 1, M'(t) = M + 1, M(0) = 0
    // M(t) = exp(t) - 1
    // M(T) = L => T = ln(L+1)
    double length = 3.0;
    LowerEnvelope env;
    env.pieces.push_back({0.0, length, Linear{1.0, 1.0}});

    auto sol = solve_ode(env, length);
    double expected_T = std::log(length + 1.0);
    EXPECT_NEAR(sol.total_T, expected_T, 1e-10);

    // Verify evaluate at midpoint
    double t_mid = expected_T / 2.0;
    double expected_M = std::exp(t_mid) - 1.0;
    EXPECT_NEAR(sol.evaluate(t_mid), expected_M, 1e-8);
}

TEST(ODESolver, SqrtQuadraticSinh) {
    // F(x) = sqrt(x^2 + 1), M(0) = 0
    // M(t) = sinh(t) (since c=0, d=1, arcsinh(0)=0)
    // M(T) = L => T = arcsinh(L)
    double length = 2.0;
    LowerEnvelope env;
    env.pieces.push_back({0.0, length, SqrtQuadratic{0.0, 1.0}});

    auto sol = solve_ode(env, length);
    double expected_T = std::asinh(length);
    EXPECT_NEAR(sol.total_T, expected_T, 1e-10);

    // Verify M(t) = sinh(t)
    double t_mid = expected_T / 2.0;
    EXPECT_NEAR(sol.evaluate(t_mid), std::sinh(t_mid), 1e-8);
}

TEST(ODESolver, MultiPieceContinuity) {
    // Two pieces: linear then sqrt-quadratic
    double x_split = 1.0;
    double length = 3.0;
    LowerEnvelope env;
    env.pieces.push_back({0.0, x_split, Linear{0.0, 1.0}});  // F=1 constant
    env.pieces.push_back({x_split, length, SqrtQuadratic{1.0, 1.0}});  // F=sqrt((x-1)^2+1)

    auto sol = solve_ode(env, length);

    // First piece: M(t)=t for t in [0,1], so piece 1 ends at t=1
    ASSERT_GE(sol.pieces.size(), 2u);
    EXPECT_NEAR(sol.pieces[0].t_end, 1.0, 1e-10);
    EXPECT_NEAR(sol.pieces[0].x_end, x_split, 1e-10);

    // Continuity: end of piece 0 == start of piece 1
    EXPECT_NEAR(sol.pieces[0].x_end, sol.pieces[1].x_start, 1e-10);
    EXPECT_NEAR(sol.pieces[0].t_end, sol.pieces[1].t_start, 1e-10);
}

TEST(ODESolver, Roundtrip) {
    double length = 5.0;
    LowerEnvelope env;
    env.pieces.push_back({0.0, length, SqrtQuadratic{2.0, 1.5}});

    auto sol = solve_ode(env, length);

    for (const auto& piece : sol.pieces) {
        EXPECT_NEAR(sol.evaluate(piece.t_start), piece.x_start, 1e-8);
        EXPECT_NEAR(sol.evaluate(piece.t_end), piece.x_end, 1e-8);
    }
}
