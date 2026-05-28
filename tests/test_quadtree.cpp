#include <gtest/gtest.h>
#include "quadtree.h"
#include <cmath>
#include <set>

using namespace af2d;

// --- BBox tests ---

TEST(BBox, ContainsPoint) {
    BBox box{0.0, 0.0, 10.0, 10.0};

    // Interior
    EXPECT_TRUE(box.contains(5.0, 5.0));

    // Boundary
    EXPECT_TRUE(box.contains(0.0, 0.0));
    EXPECT_TRUE(box.contains(10.0, 10.0));
    EXPECT_TRUE(box.contains(0.0, 5.0));
    EXPECT_TRUE(box.contains(5.0, 0.0));

    // Exterior
    EXPECT_FALSE(box.contains(-0.1, 5.0));
    EXPECT_FALSE(box.contains(5.0, -0.1));
    EXPECT_FALSE(box.contains(10.1, 5.0));
    EXPECT_FALSE(box.contains(5.0, 10.1));
}

TEST(BBox, Intersects) {
    BBox a{0.0, 0.0, 10.0, 10.0};

    // Overlapping
    BBox b{5.0, 5.0, 15.0, 15.0};
    EXPECT_TRUE(a.intersects(b));
    EXPECT_TRUE(b.intersects(a));

    // Non-overlapping
    BBox c{20.0, 20.0, 30.0, 30.0};
    EXPECT_FALSE(a.intersects(c));
    EXPECT_FALSE(c.intersects(a));

    // Touching edge
    BBox d{10.0, 0.0, 20.0, 10.0};
    EXPECT_TRUE(a.intersects(d));
    EXPECT_TRUE(d.intersects(a));

    // Touching corner
    BBox e{10.0, 10.0, 20.0, 20.0};
    EXPECT_TRUE(a.intersects(e));
}

TEST(BBox, Expanded) {
    BBox box{1.0, 2.0, 3.0, 4.0};
    BBox exp = box.expanded(0.5);
    EXPECT_DOUBLE_EQ(exp.x_min, 0.5);
    EXPECT_DOUBLE_EQ(exp.y_min, 1.5);
    EXPECT_DOUBLE_EQ(exp.x_max, 3.5);
    EXPECT_DOUBLE_EQ(exp.y_max, 4.5);
}

// --- Quadtree structural tests ---

TEST(Quadtree, EmptyTree) {
    PSLG pslg;
    Quadtree qt;
    qt.build(pslg);

    std::vector<FeatureItem> results;
    qt.query({-100, -100, 100, 100}, results);
    EXPECT_TRUE(results.empty());
}

TEST(Quadtree, SingleVertex) {
    PSLG pslg;
    pslg.vertices = {{5.0, 5.0}};

    Quadtree qt;
    qt.build(pslg);

    std::vector<FeatureItem> results;
    qt.query({4.0, 4.0, 6.0, 6.0}, results);
    ASSERT_EQ(results.size(), 1u);
    EXPECT_EQ(results[0].type, FeatureType::Vertex);
    EXPECT_EQ(results[0].index, 0);
}

TEST(Quadtree, MultipleVerticesRangeQuery) {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {10, 0}, {20, 0}, {30, 0}};

    Quadtree qt;
    qt.build(pslg);

    // Query that covers only first two vertices
    std::vector<FeatureItem> results;
    qt.query({-1, -1, 11, 1}, results);

    std::set<int> found_indices;
    for (const auto& item : results) {
        if (item.type == FeatureType::Vertex) {
            found_indices.insert(item.index);
        }
    }
    EXPECT_TRUE(found_indices.count(0));
    EXPECT_TRUE(found_indices.count(1));
    EXPECT_FALSE(found_indices.count(2));
    EXPECT_FALSE(found_indices.count(3));
}

TEST(Quadtree, SegmentBBoxQuery) {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {10, 0}};
    pslg.segments = {{0, 1}};

    Quadtree qt;
    qt.build(pslg);

    // Query overlapping the segment
    std::vector<FeatureItem> results;
    qt.query({3, -1, 7, 1}, results);

    bool found_segment = false;
    for (const auto& item : results) {
        if (item.type == FeatureType::Segment && item.index == 0) {
            found_segment = true;
        }
    }
    EXPECT_TRUE(found_segment);
}

TEST(Quadtree, QueryOutsideRange) {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    pslg.segments = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

    Quadtree qt;
    qt.build(pslg);

    std::vector<FeatureItem> results;
    qt.query({100, 100, 200, 200}, results);
    EXPECT_TRUE(results.empty());
}

TEST(Quadtree, BuildFromPSLG) {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    pslg.segments = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

    Quadtree qt;
    qt.build(pslg);

    // Query everything
    std::vector<FeatureItem> results;
    qt.query({-1, -1, 2, 2}, results);

    // Should find 4 vertices + 4 segments = 8 items
    int vertex_count = 0, segment_count = 0;
    for (const auto& item : results) {
        if (item.type == FeatureType::Vertex) vertex_count++;
        else segment_count++;
    }
    EXPECT_EQ(vertex_count, 4);
    EXPECT_EQ(segment_count, 4);
}

TEST(Quadtree, DegenerateSameLocation) {
    PSLG pslg;
    pslg.vertices = {{5.0, 5.0}, {5.0, 5.0}, {5.0, 5.0}};

    Quadtree qt;
    qt.build(pslg);

    std::vector<FeatureItem> results;
    qt.query({4, 4, 6, 6}, results);

    int vertex_count = 0;
    for (const auto& item : results) {
        if (item.type == FeatureType::Vertex) vertex_count++;
    }
    EXPECT_EQ(vertex_count, 3);
}

// --- Distance helper tests ---

TEST(DistanceHelpers, PointToSegmentPerpendicular) {
    // Point above midpoint of horizontal segment
    Point2D p{5.0, 3.0};
    Point2D a{0.0, 0.0}, b{10.0, 0.0};
    double d = min_distance_point_to_segment(p, a, b);
    EXPECT_NEAR(d, 3.0, 1e-12);
}

TEST(DistanceHelpers, PointToSegmentEndpoint) {
    // Point beyond segment end
    Point2D p{15.0, 0.0};
    Point2D a{0.0, 0.0}, b{10.0, 0.0};
    double d = min_distance_point_to_segment(p, a, b);
    EXPECT_NEAR(d, 5.0, 1e-12);
}

TEST(DistanceHelpers, PointToSegmentOnSegment) {
    Point2D p{5.0, 0.0};
    Point2D a{0.0, 0.0}, b{10.0, 0.0};
    double d = min_distance_point_to_segment(p, a, b);
    EXPECT_NEAR(d, 0.0, 1e-12);
}

TEST(DistanceHelpers, SegmentToSegmentParallel) {
    // Two parallel horizontal segments, 3 apart
    Point2D a1{0, 0}, b1{10, 0};
    Point2D a2{0, 3}, b2{10, 3};
    double d = min_distance_segment_to_segment(a1, b1, a2, b2);
    EXPECT_NEAR(d, 3.0, 1e-12);
}

TEST(DistanceHelpers, SegmentToSegmentPerpendicular) {
    // T-shape: horizontal and vertical segments
    Point2D a1{0, 0}, b1{10, 0};  // horizontal
    Point2D a2{5, 2}, b2{5, 8};   // vertical, 2 above
    double d = min_distance_segment_to_segment(a1, b1, a2, b2);
    EXPECT_NEAR(d, 2.0, 1e-12);
}

TEST(DistanceHelpers, SegmentToSegmentEndpointToEndpoint) {
    // Two collinear segments separated by gap
    Point2D a1{0, 0}, b1{3, 0};
    Point2D a2{7, 0}, b2{10, 0};
    double d = min_distance_segment_to_segment(a1, b1, a2, b2);
    EXPECT_NEAR(d, 4.0, 1e-12);
}
