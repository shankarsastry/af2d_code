#include <gtest/gtest.h>
#include "triangle_wrapper.h"
#include "segment_splitter.h"
#include <cmath>
#include <set>

using namespace af2d;

namespace {

double triangle_area(const Point2D& a, const Point2D& b, const Point2D& c) {
    return 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
}

bool edge_exists(const TriangleOutput& out, int p, int q) {
    for (const auto& tri : out.triangles) {
        for (int i = 0; i < 3; ++i) {
            int a = tri[i], b = tri[(i + 1) % 3];
            if ((a == p && b == q) || (a == q && b == p)) return true;
        }
    }
    return false;
}

} // namespace

TEST(EndToEnd, SquareFullPipeline) {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    pslg.segments = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    pslg.validate();

    auto out = triangulate(pslg);

    EXPECT_GE(out.num_triangles, 2);
    EXPECT_GE(out.num_points, 4);

    // All triangles have nonzero area
    for (const auto& tri : out.triangles) {
        double area = triangle_area(out.points[tri[0]], out.points[tri[1]], out.points[tri[2]]);
        EXPECT_GT(std::abs(area), 1e-15);
    }
}

TEST(EndToEnd, LShapeGradedMesh) {
    PSLG pslg;
    pslg.vertices = {
        {0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}
    };
    pslg.segments = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}
    };
    pslg.validate();

    auto split = split_segments(pslg);
    auto out = triangulate_pslg(split);

    // Should produce a reasonable number of triangles
    EXPECT_GE(out.num_triangles, 4);

    // All triangles oriented
    for (const auto& tri : out.triangles) {
        double area = triangle_area(out.points[tri[0]], out.points[tri[1]], out.points[tri[2]]);
        EXPECT_GT(std::abs(area), 1e-15);
    }
}

TEST(EndToEnd, NarrowChannel) {
    // Two long parallel segments close together forming a narrow channel
    PSLG pslg;
    pslg.vertices = {
        {0, 0}, {10, 0}, {10, 0.2}, {0, 0.2}
    };
    pslg.segments = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}
    };
    pslg.validate();

    auto split = split_segments(pslg);

    // The long segments should be split into many subsegments
    // to resolve the narrow width
    EXPECT_GT(split.all_segments.size(), 4u);

    auto out = triangulate_pslg(split);
    EXPECT_GE(out.num_triangles, 4);
}

TEST(EndToEnd, InputSegmentsPreserved) {
    // After triangulation, original input segments should appear as triangle edges
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    pslg.segments = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

    auto split = split_segments(pslg);
    // Request edges in output with 'e' switch
    auto out = triangulate_pslg(split, "pzQe");

    // Each sub-segment should be an edge in the triangulation
    for (const auto& seg : split.all_segments) {
        bool found = edge_exists(out, seg.p, seg.q);
        EXPECT_TRUE(found) << "Segment (" << seg.p << ", " << seg.q
                           << ") not found as triangle edge";
    }
}
