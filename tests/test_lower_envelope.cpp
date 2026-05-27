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
