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

TEST(ODESolver, EvaluateAtBoundaries) {
    // Multi-piece solution: evaluate at exact t_start and t_end of each piece
    double length = 6.0;
    LowerEnvelope env;
    env.pieces.push_back({0.0, 2.0, Linear{0.0, 1.0}});       // F=1 constant
    env.pieces.push_back({2.0, 4.0, SqrtQuadratic{3.0, 0.5}}); // sqrt((x-3)^2+0.25)
    env.pieces.push_back({4.0, 6.0, Linear{0.5, 2.0}});        // F=0.5x+2

    auto sol = solve_ode(env, length);
    ASSERT_GE(sol.pieces.size(), 2u);

    for (const auto& piece : sol.pieces) {
        EXPECT_NEAR(sol.evaluate(piece.t_start), piece.x_start, 1e-8)
            << "Mismatch at t_start=" << piece.t_start;
        EXPECT_NEAR(sol.evaluate(piece.t_end), piece.x_end, 1e-8)
            << "Mismatch at t_end=" << piece.t_end;
    }
}

TEST(ODESolver, ManyPieces) {
    // Construct envelope with 10+ linear pieces, solve ODE, verify evaluate
    double length = 10.0;
    LowerEnvelope env;
    int n_pieces = 12;
    double dx = length / n_pieces;
    for (int i = 0; i < n_pieces; ++i) {
        double x0 = i * dx;
        double x1 = (i + 1) * dx;
        // Alternating constant values: 1.0 and 2.0
        double c = (i % 2 == 0) ? 1.0 : 2.0;
        env.pieces.push_back({x0, x1, Linear{0.0, c}});
    }

    auto sol = solve_ode(env, length);
    ASSERT_GE(sol.pieces.size(), 10u);
    EXPECT_NEAR(sol.evaluate(0.0), 0.0, 1e-10);
    EXPECT_NEAR(sol.evaluate(sol.total_T), length, 1e-8);

    // Evaluate at several intermediate t values
    for (int k = 1; k < 20; ++k) {
        double t = k * sol.total_T / 20.0;
        double x = sol.evaluate(t);
        EXPECT_GE(x, 0.0 - 1e-10);
        EXPECT_LE(x, length + 1e-10);
    }

    // Verify monotonicity
    double prev_x = 0.0;
    for (int k = 1; k <= 100; ++k) {
        double t = k * sol.total_T / 100.0;
        double x = sol.evaluate(t);
        EXPECT_GE(x, prev_x - 1e-10) << "Non-monotonic at t=" << t;
        prev_x = x;
    }
}
