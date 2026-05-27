#include <gtest/gtest.h>
#include "triangle_wrapper.h"
#include <cmath>

using namespace af2d;

namespace {

double triangle_area(const Point2D& a, const Point2D& b, const Point2D& c) {
    return 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
}

} // namespace

TEST(TriangleWrapper, SquareBasic) {
    SplitResult sr;
    sr.all_vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    sr.all_segments = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

    auto out = triangulate_pslg(sr);
    EXPECT_EQ(out.num_points, 4);
    EXPECT_EQ(out.num_triangles, 2);
}

TEST(TriangleWrapper, TrianglesPositiveArea) {
    SplitResult sr;
    sr.all_vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    sr.all_segments = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

    auto out = triangulate_pslg(sr);
    for (const auto& tri : out.triangles) {
        double area = triangle_area(out.points[tri[0]], out.points[tri[1]], out.points[tri[2]]);
        EXPECT_GT(std::abs(area), 1e-15) << "Degenerate triangle found";
    }
}

TEST(TriangleWrapper, PreSplitSquare) {
    // Split a square first, then triangulate
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    pslg.segments = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

    auto split = split_segments(pslg);
    auto out = triangulate_pslg(split);

    // Should produce a valid triangulation
    EXPECT_GE(out.num_points, 4);
    EXPECT_GE(out.num_triangles, 2);

    // All triangles should have positive area
    for (const auto& tri : out.triangles) {
        double area = triangle_area(out.points[tri[0]], out.points[tri[1]], out.points[tri[2]]);
        EXPECT_GT(std::abs(area), 1e-15);
    }
}

TEST(TriangleWrapper, ConvenienceFunction) {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    pslg.segments = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

    auto out = triangulate(pslg);
    EXPECT_GE(out.num_triangles, 2);
}

TEST(TriangleWrapper, SquareWithHole) {
    // Outer square with inner square hole
    SplitResult sr;
    sr.all_vertices = {
        // Outer square
        {0, 0}, {4, 0}, {4, 4}, {0, 4},
        // Inner square (hole boundary)
        {1, 1}, {3, 1}, {3, 3}, {1, 3}
    };
    sr.all_segments = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},  // outer
        {4, 5}, {5, 6}, {6, 7}, {7, 4}   // inner (hole)
    };

    // Use "pzQ" with a hole point
    // For holes, we need to use Triangle's hole mechanism
    // For now, just verify it triangulates without the hole spec
    auto out = triangulate_pslg(sr, "pzQ");
    EXPECT_GE(out.num_triangles, 4);
}
