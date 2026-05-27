#include <gtest/gtest.h>
#include "distance.h"
#include <cmath>

using namespace af2d;

TEST(Distance, VertexAboveMidpoint) {
    // Segment along x-axis from (0,0) to (2,0), vertex at (1,1)
    PSLG pslg;
    pslg.vertices = {{0, 0}, {2, 0}, {1, 1}};
    pslg.segments = {{0, 1}};

    auto df = distance_to_vertex(pslg, 0, 2);
    ASSERT_EQ(df.pieces.size(), 1u);

    // At x=1 (directly below vertex): distance = 1
    EXPECT_NEAR(df.evaluate(1.0), 1.0, 1e-12);

    // At x=0: distance = sqrt(1+1) = sqrt(2)
    EXPECT_NEAR(df.evaluate(0.0), std::sqrt(2.0), 1e-12);

    // At x=2: distance = sqrt(1+1) = sqrt(2)
    EXPECT_NEAR(df.evaluate(2.0), std::sqrt(2.0), 1e-12);
}

TEST(Distance, VertexOnExtension) {
    // Segment (0,0) to (1,0), vertex at (2,0) on the extension
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}, {2, 0}};
    pslg.segments = {{0, 1}};

    auto df = distance_to_vertex(pslg, 0, 2);

    // At x=1 (end of segment): distance to (2,0) = 1
    EXPECT_NEAR(df.evaluate(1.0), 1.0, 1e-12);

    // At x=0: distance = 2
    EXPECT_NEAR(df.evaluate(0.0), 2.0, 1e-12);
}

TEST(Distance, ParallelSegments) {
    // S_i: (0,0) to (2,0), S_j: (0,1) to (2,1)
    PSLG pslg;
    pslg.vertices = {{0, 0}, {2, 0}, {0, 1}, {2, 1}};
    pslg.segments = {{0, 1}, {2, 3}};

    auto df = distance_to_segment(pslg, 0, 1);

    // Parallel segments at distance 1 -> constant perpendicular distance
    EXPECT_NEAR(df.evaluate(0.0), 1.0, 1e-6);
    EXPECT_NEAR(df.evaluate(1.0), 1.0, 1e-6);
    EXPECT_NEAR(df.evaluate(2.0), 1.0, 1e-6);
}

TEST(Distance, PerpendicularSegments) {
    // S_i: (0,0) to (4,0), S_j: (2,1) to (2,3) (perpendicular, above midpoint)
    PSLG pslg;
    pslg.vertices = {{0, 0}, {4, 0}, {2, 1}, {2, 3}};
    pslg.segments = {{0, 1}, {2, 3}};

    auto df = distance_to_segment(pslg, 0, 1);

    // At x=2: closest point on S_j is (2,1), distance = 1
    EXPECT_NEAR(df.evaluate(2.0), 1.0, 1e-6);

    // At x=0: closest point on S_j is (2,1), distance = sqrt(4+1) = sqrt(5)
    EXPECT_NEAR(df.evaluate(0.0), std::sqrt(5.0), 1e-6);
}

TEST(Distance, EndpointProximity) {
    // S_i: (0,0) to (1,0), S_j: (1.5,0.5) to (2.5,0.5)
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}, {1.5, 0.5}, {2.5, 0.5}};
    pslg.segments = {{0, 1}, {2, 3}};

    auto df = distance_to_segment(pslg, 0, 1);

    // At x=1 (nearest end): distance should be to nearest point on S_j = (1.5,0.5)
    // distance = sqrt(0.25 + 0.25) = sqrt(0.5)
    EXPECT_NEAR(df.evaluate(1.0), std::sqrt(0.5), 1e-6);
}

TEST(Distance, SqrtQuadraticEvaluate) {
    SqrtQuadratic sq{3.0, 4.0};
    // sqrt((5-3)^2 + 16) = sqrt(4+16) = sqrt(20)
    EXPECT_NEAR(sq.evaluate(5.0), std::sqrt(20.0), 1e-12);
}

TEST(Distance, LinearEvaluate) {
    Linear lin{2.0, 3.0};
    EXPECT_NEAR(lin.evaluate(1.0), 5.0, 1e-12);
}

TEST(Distance, FarthestEndpointBasic) {
    // Segment of length 4: max(x, 4-x)
    auto df = distance_to_farthest_endpoint(4.0);
    ASSERT_EQ(df.pieces.size(), 2u);

    // At x=0: max(0, 4) = 4
    EXPECT_NEAR(df.evaluate(0.0), 4.0, 1e-12);

    // At x=2 (midpoint): max(2, 2) = 2
    EXPECT_NEAR(df.evaluate(2.0), 2.0, 1e-12);

    // At x=4: max(4, 0) = 4
    EXPECT_NEAR(df.evaluate(4.0), 4.0, 1e-12);

    // At x=1: max(1, 3) = 3
    EXPECT_NEAR(df.evaluate(1.0), 3.0, 1e-12);

    // At x=3: max(3, 1) = 3
    EXPECT_NEAR(df.evaluate(3.0), 3.0, 1e-12);
}

TEST(Distance, FarthestEndpointUnitLength) {
    auto df = distance_to_farthest_endpoint(1.0);
    ASSERT_EQ(df.pieces.size(), 2u);

    // At endpoints: value = 1
    EXPECT_NEAR(df.evaluate(0.0), 1.0, 1e-12);
    EXPECT_NEAR(df.evaluate(1.0), 1.0, 1e-12);

    // At midpoint: value = 0.5
    EXPECT_NEAR(df.evaluate(0.5), 0.5, 1e-12);
}

TEST(Distance, FarthestEndpointCapsLFS) {
    // For a unit square, vertex (1,1) is at distance sqrt((x-1)^2+1) from
    // the bottom edge (0,0)-(1,0). At x=0 this is sqrt(2) ≈ 1.414.
    // The farthest-endpoint function at x=0 is max(0, 1) = 1, which is smaller.
    // So the farthest-endpoint function should dominate near endpoints.
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}, {1, 1}};
    pslg.segments = {{0, 1}};

    auto df_vertex = distance_to_vertex(pslg, 0, 2);
    auto df_farthest = distance_to_farthest_endpoint(1.0);

    // At x=0: vertex distance = sqrt(2), farthest endpoint = 1
    EXPECT_GT(df_vertex.evaluate(0.0), df_farthest.evaluate(0.0));

    // At x=0.5: vertex distance = sqrt(0.25+1) ≈ 1.118, farthest endpoint = 0.5
    EXPECT_GT(df_vertex.evaluate(0.5), df_farthest.evaluate(0.5));
}
