#include <gtest/gtest.h>
#include "io.h"
#include <cmath>
#include <filesystem>
#include <fstream>

using namespace af2d;

namespace fs = std::filesystem;

// Path to test data relative to the test executable's working directory.
// CMake/CTest runs from the build dir, so we need the source-tree path.
static const std::string TEST_DATA = CMAKE_SOURCE_DIR "/tests/test_data";

// --- Read tests ---

TEST(IO, ReadNodeBasic) {
    auto pslg = read_node(TEST_DATA + "/square.node");
    ASSERT_EQ(pslg.vertices.size(), 4u);
    EXPECT_DOUBLE_EQ(pslg.vertices[0].x, 0.0);
    EXPECT_DOUBLE_EQ(pslg.vertices[0].y, 0.0);
    EXPECT_DOUBLE_EQ(pslg.vertices[1].x, 1.0);
    EXPECT_DOUBLE_EQ(pslg.vertices[1].y, 0.0);
    EXPECT_DOUBLE_EQ(pslg.vertices[2].x, 1.0);
    EXPECT_DOUBLE_EQ(pslg.vertices[2].y, 1.0);
    EXPECT_DOUBLE_EQ(pslg.vertices[3].x, 0.0);
    EXPECT_DOUBLE_EQ(pslg.vertices[3].y, 1.0);
    EXPECT_TRUE(pslg.segments.empty());
    EXPECT_TRUE(pslg.holes.empty());
}

TEST(IO, ReadPolyBasic) {
    auto pslg = read_poly(TEST_DATA + "/square.poly");
    ASSERT_EQ(pslg.vertices.size(), 4u);
    ASSERT_EQ(pslg.segments.size(), 4u);
    EXPECT_EQ(pslg.segments[0].p, 0);
    EXPECT_EQ(pslg.segments[0].q, 1);
    EXPECT_EQ(pslg.segments[3].p, 3);
    EXPECT_EQ(pslg.segments[3].q, 0);
    EXPECT_TRUE(pslg.holes.empty());
}

TEST(IO, ReadPolyWithHole) {
    auto pslg = read_poly(TEST_DATA + "/square_with_hole.poly");
    ASSERT_EQ(pslg.vertices.size(), 8u);
    ASSERT_EQ(pslg.segments.size(), 8u);
    ASSERT_EQ(pslg.holes.size(), 1u);
    EXPECT_DOUBLE_EQ(pslg.holes[0].x, 2.0);
    EXPECT_DOUBLE_EQ(pslg.holes[0].y, 2.0);
}

TEST(IO, ReadPolyOneIndexed) {
    auto pslg = read_poly(TEST_DATA + "/1indexed.poly");
    ASSERT_EQ(pslg.vertices.size(), 4u);
    ASSERT_EQ(pslg.segments.size(), 4u);
    // Segments should be remapped to 0-based
    EXPECT_EQ(pslg.segments[0].p, 0);
    EXPECT_EQ(pslg.segments[0].q, 1);
    EXPECT_EQ(pslg.segments[3].p, 3);
    EXPECT_EQ(pslg.segments[3].q, 0);
}

TEST(IO, ReadPolyCommentsAndBlanks) {
    auto pslg = read_poly(TEST_DATA + "/comments.poly");
    ASSERT_EQ(pslg.vertices.size(), 4u);
    ASSERT_EQ(pslg.segments.size(), 4u);
    EXPECT_TRUE(pslg.holes.empty());
}

TEST(IO, ReadPolyMissingFile) {
    EXPECT_THROW(read_poly(TEST_DATA + "/nonexistent.poly"), std::runtime_error);
}

TEST(IO, ReadNodeFromPoly) {
    // .poly with 0 vertices should read from companion .node file
    auto pslg = read_poly(TEST_DATA + "/zero_vertices.poly");
    ASSERT_EQ(pslg.vertices.size(), 4u);
    ASSERT_EQ(pslg.segments.size(), 4u);
    EXPECT_DOUBLE_EQ(pslg.vertices[0].x, 0.0);
    EXPECT_DOUBLE_EQ(pslg.vertices[1].x, 1.0);
}

TEST(IO, ReadNodeFromPolyOneIndexed) {
    // .poly with 0 vertices and 1-based segment indices + companion .node
    auto pslg = read_poly(TEST_DATA + "/1indexed_companion.poly");
    ASSERT_EQ(pslg.vertices.size(), 4u);
    ASSERT_EQ(pslg.segments.size(), 4u);
    // Segments should be remapped to 0-based
    EXPECT_EQ(pslg.segments[0].p, 0);
    EXPECT_EQ(pslg.segments[0].q, 1);
    EXPECT_EQ(pslg.segments[3].p, 3);
    EXPECT_EQ(pslg.segments[3].q, 0);
}

// --- Write tests ---

TEST(IO, WriteNodeRoundtrip) {
    TriangleOutput out;
    out.points = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    out.num_points = 4;

    auto tmp = fs::temp_directory_path() / "af2d_test_write.node";
    write_node(tmp.string(), out);

    auto pslg = read_node(tmp.string());
    ASSERT_EQ(pslg.vertices.size(), 4u);
    EXPECT_DOUBLE_EQ(pslg.vertices[0].x, 0.0);
    EXPECT_DOUBLE_EQ(pslg.vertices[1].x, 1.0);
    EXPECT_DOUBLE_EQ(pslg.vertices[2].y, 1.0);

    fs::remove(tmp);
}

TEST(IO, WriteEleBasic) {
    TriangleOutput out;
    out.triangles = {{0, 1, 2}, {0, 2, 3}};
    out.num_triangles = 2;

    auto tmp = fs::temp_directory_path() / "af2d_test_write.ele";
    write_ele(tmp.string(), out);

    // Read back and verify format
    std::ifstream f(tmp.string());
    int nt, corners, attr;
    f >> nt >> corners >> attr;
    EXPECT_EQ(nt, 2);
    EXPECT_EQ(corners, 3);

    int id, v0, v1, v2;
    f >> id >> v0 >> v1 >> v2;
    EXPECT_EQ(id, 0);
    EXPECT_EQ(v0, 0);
    EXPECT_EQ(v1, 1);
    EXPECT_EQ(v2, 2);

    fs::remove(tmp);
}

TEST(IO, WriteEdgeBasic) {
    TriangleOutput out;
    out.edges = {{0, 1}, {1, 2}, {2, 0}};
    out.num_edges = 3;

    auto tmp = fs::temp_directory_path() / "af2d_test_write.edge";
    write_edge(tmp.string(), out);

    std::ifstream f(tmp.string());
    int ne, bm;
    f >> ne >> bm;
    EXPECT_EQ(ne, 3);
    EXPECT_EQ(bm, 0);

    int id, p, q;
    f >> id >> p >> q;
    EXPECT_EQ(id, 0);
    EXPECT_EQ(p, 0);
    EXPECT_EQ(q, 1);

    fs::remove(tmp);
}

// --- Integration test ---

TEST(IO, ReadSplitTriangulateWrite) {
    auto pslg = read_poly(TEST_DATA + "/square.poly");
    pslg.validate();

    auto split = split_segments(pslg);
    auto out = triangulate_pslg(split, "pzQe");

    ASSERT_GE(out.num_points, 4);
    ASSERT_GE(out.num_triangles, 2);
    ASSERT_GE(out.num_edges, 1);

    auto tmp_dir = fs::temp_directory_path() / "af2d_test_integration";
    fs::create_directories(tmp_dir);
    std::string base = (tmp_dir / "square.1").string();

    write_node(base + ".node", out);
    write_ele(base + ".ele", out);
    write_edge(base + ".edge", out);

    // Read back .node and verify point count
    auto readback = read_node(base + ".node");
    EXPECT_EQ(static_cast<int>(readback.vertices.size()), out.num_points);

    fs::remove_all(tmp_dir);
}

// --- Holes integration test ---

TEST(TriangleWrapper, HolesPassedThrough) {
    // Test that holes propagate through SplitResult to Triangle.
    // Use a SplitResult directly (no splitting) to isolate the hole-passing logic.
    SplitResult sr;
    sr.all_vertices = {
        {0, 0}, {4, 0}, {4, 4}, {0, 4},
        {1, 1}, {3, 1}, {3, 3}, {1, 3}
    };
    sr.all_segments = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4}
    };
    sr.holes = {{2.0, 2.0}};

    auto out = triangulate_pslg(sr, "pzQe");
    EXPECT_GE(out.num_triangles, 4);

    // With hole, total area should be approximately 16 - 4 = 12
    double total_area = 0.0;
    for (const auto& tri : out.triangles) {
        const auto& a = out.points[tri[0]];
        const auto& b = out.points[tri[1]];
        const auto& c = out.points[tri[2]];
        total_area += std::abs(0.5 * ((b.x - a.x) * (c.y - a.y) -
                                       (c.x - a.x) * (b.y - a.y)));
    }
    EXPECT_NEAR(total_area, 12.0, 0.1);
}

TEST(TriangleWrapper, HolesPropagatedThroughSplit) {
    // Test that holes flow through split_segments into SplitResult.
    PSLG pslg;
    pslg.vertices = {
        {0, 0}, {4, 0}, {4, 4}, {0, 4},
        {1, 1}, {3, 1}, {3, 3}, {1, 3}
    };
    pslg.segments = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4}
    };
    pslg.holes = {{2.0, 2.0}};

    auto split = split_segments(pslg);
    ASSERT_EQ(split.holes.size(), 1u);
    EXPECT_DOUBLE_EQ(split.holes[0].x, 2.0);
    EXPECT_DOUBLE_EQ(split.holes[0].y, 2.0);
}
