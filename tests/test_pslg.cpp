#include <gtest/gtest.h>
#include "pslg.h"

using namespace af2d;

namespace {

PSLG make_unit_square() {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    pslg.segments = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    return pslg;
}

} // namespace

TEST(PSLG, SquareVertexCount) {
    auto pslg = make_unit_square();
    EXPECT_EQ(pslg.vertices.size(), 4u);
}

TEST(PSLG, SquareSegmentCount) {
    auto pslg = make_unit_square();
    EXPECT_EQ(pslg.segments.size(), 4u);
}

TEST(PSLG, SegmentLengths) {
    auto pslg = make_unit_square();
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(pslg.segment_length(i), 1.0, 1e-12);
    }
}

TEST(PSLG, PointOnSegmentEndpoints) {
    auto pslg = make_unit_square();
    // Segment 0: (0,0) -> (1,0)
    auto p0 = pslg.point_on_segment(0, 0.0);
    EXPECT_NEAR(p0.x, 0.0, 1e-12);
    EXPECT_NEAR(p0.y, 0.0, 1e-12);

    auto p1 = pslg.point_on_segment(0, 1.0);
    EXPECT_NEAR(p1.x, 1.0, 1e-12);
    EXPECT_NEAR(p1.y, 0.0, 1e-12);
}

TEST(PSLG, PointOnSegmentMidpoint) {
    auto pslg = make_unit_square();
    auto pm = pslg.point_on_segment(0, 0.5);
    EXPECT_NEAR(pm.x, 0.5, 1e-12);
    EXPECT_NEAR(pm.y, 0.0, 1e-12);
}

TEST(PSLG, AdjacencySegments) {
    auto pslg = make_unit_square();
    // Segments 0 and 1 share vertex 1
    EXPECT_TRUE(pslg.is_adjacent(0, 1));
    // Segments 0 and 2 share no vertices
    EXPECT_FALSE(pslg.is_adjacent(0, 2));
    // Segments 0 and 3 share vertex 0
    EXPECT_TRUE(pslg.is_adjacent(0, 3));
}

TEST(PSLG, AdjacencyVertex) {
    auto pslg = make_unit_square();
    EXPECT_TRUE(pslg.is_adjacent_vertex(0, 0));
    EXPECT_TRUE(pslg.is_adjacent_vertex(0, 1));
    EXPECT_FALSE(pslg.is_adjacent_vertex(0, 2));
    EXPECT_FALSE(pslg.is_adjacent_vertex(0, 3));
}

TEST(PSLG, ValidateAcceptsGood) {
    auto pslg = make_unit_square();
    EXPECT_NO_THROW(pslg.validate());
}

TEST(PSLG, ValidateRejectsOutOfRange) {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}};
    pslg.segments = {{0, 5}};
    EXPECT_THROW(pslg.validate(), std::invalid_argument);
}

TEST(PSLG, ValidateRejectsSameEndpoints) {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}};
    pslg.segments = {{0, 0}};
    EXPECT_THROW(pslg.validate(), std::invalid_argument);
}

TEST(PSLG, ValidateRejectsZeroLength) {
    PSLG pslg;
    pslg.vertices = {{1, 2}, {1, 2}};
    pslg.segments = {{0, 1}};
    EXPECT_THROW(pslg.validate(), std::invalid_argument);
}
