#include <gtest/gtest.h>
#include "lower_envelope.h"
#include <cmath>

using namespace af2d;

TEST(LowerEnvelope, SingleFunction) {
    DistanceFunction df;
    df.pieces.push_back({0.0, 2.0, SqrtQuadratic{1.0, 1.0}});

    auto env = compute_lower_envelope({df}, 2.0);
    ASSERT_EQ(env.pieces.size(), 1u);
    EXPECT_NEAR(env.evaluate(1.0), 1.0, 1e-10);  // sqrt((1-1)^2 + 1) = 1
}

TEST(LowerEnvelope, TwoSymmetricVertices) {
    // Two vertices equidistant from segment endpoints
    // f1(x) = sqrt((x-0)^2 + 1) = sqrt(x^2+1), vertex at (0,1) above left end
    // f2(x) = sqrt((x-2)^2 + 1), vertex at (2,1) above right end
    // They cross at x=1 (by symmetry)
    DistanceFunction df1, df2;
    df1.pieces.push_back({0.0, 2.0, SqrtQuadratic{0.0, 1.0}});
    df2.pieces.push_back({0.0, 2.0, SqrtQuadratic{2.0, 1.0}});

    auto env = compute_lower_envelope({df1, df2}, 2.0);

    // At x=0: f1=1, f2=sqrt(5) -> f1 wins
    EXPECT_NEAR(env.evaluate(0.0), 1.0, 1e-6);
    // At x=2: f1=sqrt(5), f2=1 -> f2 wins
    EXPECT_NEAR(env.evaluate(2.0), 1.0, 1e-6);
    // At x=1: both = sqrt(2), should be sqrt(2)
    EXPECT_NEAR(env.evaluate(1.0), std::sqrt(2.0), 1e-6);
    // Should have 2 pieces (switching at x=1)
    EXPECT_EQ(env.pieces.size(), 2u);
}

TEST(LowerEnvelope, OneDominantFunction) {
    // f1(x) = sqrt(x^2+1), f2(x) = sqrt(x^2+100) -> f1 always less
    DistanceFunction df1, df2;
    df1.pieces.push_back({0.0, 2.0, SqrtQuadratic{0.0, 1.0}});
    df2.pieces.push_back({0.0, 2.0, SqrtQuadratic{0.0, 10.0}});

    auto env = compute_lower_envelope({df1, df2}, 2.0);
    EXPECT_EQ(env.pieces.size(), 1u);
}

TEST(LowerEnvelope, ThreeVertices) {
    // f1 at (0,1), f2 at (1,0.5), f3 at (2,1)
    DistanceFunction df1, df2, df3;
    df1.pieces.push_back({0.0, 2.0, SqrtQuadratic{0.0, 1.0}});
    df2.pieces.push_back({0.0, 2.0, SqrtQuadratic{1.0, 0.5}});
    df3.pieces.push_back({0.0, 2.0, SqrtQuadratic{2.0, 1.0}});

    auto env = compute_lower_envelope({df1, df2, df3}, 2.0);

    // f2 should dominate near x=1
    EXPECT_NEAR(env.evaluate(1.0), 0.5, 1e-6);
    // At x=0: f1=1, f2=sqrt(1+0.25)=1.118, f3=sqrt(5) -> f1 wins
    EXPECT_NEAR(env.evaluate(0.0), 1.0, 1e-6);
}

TEST(LowerEnvelope, LinearSqrtQuadraticCrossing) {
    // f1(x) = 0.5*x + 0.5 (linear), f2(x) = sqrt((x-1)^2+0.01) (sqrt-quad, close vertex)
    DistanceFunction df1, df2;
    df1.pieces.push_back({0.0, 2.0, Linear{0.5, 0.5}});
    df2.pieces.push_back({0.0, 2.0, SqrtQuadratic{1.0, 0.1}});

    auto env = compute_lower_envelope({df1, df2}, 2.0);

    // At x=1: f1=1.0, f2=0.1 -> f2 wins
    EXPECT_NEAR(env.evaluate(1.0), 0.1, 1e-6);
    // At x=0: f1=0.5, f2=sqrt(1+0.01)≈1.005 -> f1 wins
    EXPECT_NEAR(env.evaluate(0.0), 0.5, 1e-6);
}

TEST(LowerEnvelope, EmptyInput) {
    auto env = compute_lower_envelope({}, 2.0);
    EXPECT_TRUE(env.pieces.empty());
}

TEST(LowerEnvelope, ManyFunctions) {
    // 12 sqrt-quadratic functions with vertices spaced along a line above the segment
    double length = 10.0;
    std::vector<DistanceFunction> funcs;
    for (int i = 0; i < 12; ++i) {
        DistanceFunction df;
        double c = i * (length / 11.0);
        double d = 0.5 + 0.1 * i;  // varying heights
        df.pieces.push_back({0.0, length, SqrtQuadratic{c, d}});
        funcs.push_back(df);
    }

    auto env = compute_lower_envelope(funcs, length);

    // Envelope should have at least 2 pieces (functions cross)
    EXPECT_GE(env.pieces.size(), 2u);
    // Envelope should not exceed number of functions
    EXPECT_LE(env.pieces.size(), 24u);

    // Verify envelope is at or below every function at sample points
    for (double x = 0.0; x <= length; x += 0.5) {
        double env_val = env.evaluate(x);
        for (const auto& f : funcs) {
            EXPECT_LE(env_val, f.evaluate(x) + 1e-6)
                << "Envelope exceeds function at x=" << x;
        }
    }
}

TEST(LowerEnvelope, AllIdenticalFunctions) {
    // 5 identical functions should produce 1 merged piece
    std::vector<DistanceFunction> funcs;
    for (int i = 0; i < 5; ++i) {
        DistanceFunction df;
        df.pieces.push_back({0.0, 3.0, SqrtQuadratic{1.5, 1.0}});
        funcs.push_back(df);
    }

    auto env = compute_lower_envelope(funcs, 3.0);
    EXPECT_EQ(env.pieces.size(), 1u);
    EXPECT_NEAR(env.evaluate(1.5), 1.0, 1e-10);
}

TEST(LowerEnvelope, MergeAdjacentPieces) {
    // Two functions where one dominates everywhere except a tiny region.
    // f1(x) = sqrt(x^2 + 1) — minimum at x=0
    // f2(x) = sqrt((x-4)^2 + 1) — minimum at x=4
    // On [0,4], they cross at x=2. Both pieces from f1 (left) and f2 (right)
    // should be single merged pieces, not fragmented.
    DistanceFunction df1, df2;
    df1.pieces.push_back({0.0, 4.0, SqrtQuadratic{0.0, 1.0}});
    df2.pieces.push_back({0.0, 4.0, SqrtQuadratic{4.0, 1.0}});

    auto env = compute_lower_envelope({df1, df2}, 4.0);
    EXPECT_EQ(env.pieces.size(), 2u);
}
