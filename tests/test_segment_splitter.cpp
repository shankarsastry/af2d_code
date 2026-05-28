#include <gtest/gtest.h>
#include "segment_splitter.h"
#include <cmath>
#include <set>

using namespace af2d;

namespace {

PSLG make_unit_square() {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    pslg.segments = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    return pslg;
}

bool point_on_line_segment(Point2D p, Point2D a, Point2D b, double tol = 1e-8) {
    double len = a.distance(b);
    if (len < 1e-15) return p.distance(a) < tol;
    double t = (p - a).dot(b - a) / (len * len);
    if (t < -tol || t > 1.0 + tol) return false;
    Point2D proj = a + (b - a) * t;
    return p.distance(proj) < tol;
}

} // namespace

// --- compute_n_star_from_angle tests ---

TEST(NStarFromAngle, Below30ReturnsAtLeastTwo) {
    for (double deg : {5.0, 10.0, 15.0, 20.0, 25.0}) {
        auto result = compute_n_star_from_angle(deg);
        EXPECT_GE(result.n_star, 2) << "at " << deg << " degrees";
        EXPECT_TRUE(result.message.empty()) << "at " << deg << " degrees";
    }
}

TEST(NStarFromAngle, TwentyDegreesIsReasonable) {
    auto result = compute_n_star_from_angle(20.0);
    EXPECT_GE(result.n_star, 2);
    EXPECT_LE(result.n_star, 30);
    EXPECT_TRUE(result.message.empty());
}

TEST(NStarFromAngle, ThirtyDegreesReturnsOneWithMessage) {
    auto result = compute_n_star_from_angle(30.0);
    EXPECT_EQ(result.n_star, 1);
    EXPECT_FALSE(result.message.empty());
}

TEST(NStarFromAngle, AboveThirtyReturnsOneWithMessage) {
    for (double deg : {31.0, 33.0, 45.0, 60.0}) {
        auto result = compute_n_star_from_angle(deg);
        EXPECT_EQ(result.n_star, 1) << "at " << deg << " degrees";
        EXPECT_FALSE(result.message.empty()) << "at " << deg << " degrees";
    }
}

TEST(NStarFromAngle, InvalidAngleThrows) {
    EXPECT_THROW(compute_n_star_from_angle(0.0), std::invalid_argument);
    EXPECT_THROW(compute_n_star_from_angle(90.0), std::invalid_argument);
    EXPECT_THROW(compute_n_star_from_angle(-5.0), std::invalid_argument);
}

// --- Default behavior: n_star=1 ---

TEST(SegmentSplitter, DefaultNStarIsOne) {
    // With n_star=1, the shortest reference segment gets 1 subsegment (no split).
    // Other segments get floor(1 * T_i / T_min) which may be > 1 if T_i > T_min.
    auto pslg = make_unit_square();
    auto result = split_segments(pslg);

    // Unit square: all segments have the same T_i = T_min, so all get 1 subsegment.
    EXPECT_EQ(result.all_segments.size(), 4u);
    EXPECT_EQ(result.all_vertices.size(), 4u);
}

TEST(SegmentSplitter, DefaultNStarProportionalSplitting) {
    // Even with n_star=1, a segment with T_i >> T_min should get more subsegments.
    PSLG pslg;
    pslg.vertices = {{0, 0}, {10, 0}, {0, 0.1}, {10, 0.1}};
    pslg.segments = {{0, 1}, {2, 3}};

    auto result = split_segments(pslg);
    // Both segments are the same length and have the same T, so n_i = 1 for each.
    // With n_star=1, T_i/T_min = 1, so each gets 1 subsegment.
    EXPECT_EQ(result.all_segments.size(), 2u);
}

// --- Explicit n_star tests ---

TEST(SegmentSplitter, ExplicitNStarFour) {
    auto pslg = make_unit_square();

    SplitParams params;
    params.n_star = 4;
    auto result = split_segments(pslg, params);

    // All 4 segments are symmetric (same T_i = T_min), each gets 4 subsegments.
    EXPECT_EQ(result.all_segments.size(), 16u);
}

TEST(SegmentSplitter, ExplicitNStarProportionalScaling) {
    // Two parallel segments of different lengths.
    PSLG pslg;
    pslg.vertices = {{0, 0}, {2, 0}, {0, 1}, {10, 1}};
    pslg.segments = {{0, 1}, {2, 3}};

    SplitParams params;
    params.n_star = 3;
    auto result = split_segments(pslg, params);

    int count_seg0 = 0, count_seg1 = 0;
    for (const auto& seg : result.all_segments) {
        const auto& a = result.all_vertices[seg.p];
        const auto& b = result.all_vertices[seg.q];
        Point2D mid = (a + b) * 0.5;
        if (mid.y < 0.5) count_seg0++;
        else count_seg1++;
    }

    // Shorter segment gets n_star=3 subsegments
    EXPECT_GE(count_seg0, 3);
    // Longer segment gets proportionally more
    EXPECT_GT(count_seg1, count_seg0);
}

TEST(SegmentSplitter, IncreasingNStarMoreSegments) {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {10, 0}, {0, 0.5}, {10, 0.5}};
    pslg.segments = {{0, 1}, {2, 3}};

    SplitParams p2, p5, p10;
    p2.n_star = 2;
    p5.n_star = 5;
    p10.n_star = 10;

    auto r2 = split_segments(pslg, p2);
    auto r5 = split_segments(pslg, p5);
    auto r10 = split_segments(pslg, p10);

    EXPECT_LE(r2.all_segments.size(), r5.all_segments.size());
    EXPECT_LE(r5.all_segments.size(), r10.all_segments.size());
}

TEST(SegmentSplitter, IsolatedSegmentSplitsByNStar) {
    // An isolated segment still has features (its own endpoints define lfs).
    // With only one segment, T_i = T_min, so n_i = n_star.
    PSLG pslg;
    pslg.vertices = {{0, 0}, {5, 0}};
    pslg.segments = {{0, 1}};

    SplitParams params;
    params.n_star = 10;
    auto result = split_segments(pslg, params);

    EXPECT_EQ(result.all_segments.size(), 10u);
    EXPECT_EQ(result.all_vertices.size(), 11u);
}

TEST(SegmentSplitter, SplitPointsOnParent) {
    auto pslg = make_unit_square();

    SplitParams params;
    params.n_star = 5;
    auto result = split_segments(pslg, params);

    int orig_count = static_cast<int>(pslg.vertices.size());
    for (int i = orig_count; i < static_cast<int>(result.all_vertices.size()); ++i) {
        const auto& pt = result.all_vertices[i];
        bool on_some_segment = false;
        for (int s = 0; s < static_cast<int>(pslg.segments.size()); ++s) {
            const auto& a = pslg.vertices[pslg.segments[s].p];
            const auto& b = pslg.vertices[pslg.segments[s].q];
            if (point_on_line_segment(pt, a, b)) {
                on_some_segment = true;
                break;
            }
        }
        EXPECT_TRUE(on_some_segment) << "Split point (" << pt.x << ", " << pt.y
                                     << ") not on any original segment";
    }
}

// --- Angle-derived n_star tests ---

TEST(SegmentSplitter, AngleDerivedNStar) {
    // Compute n_star from angle, then use it.
    auto pslg = make_unit_square();

    auto nsr = compute_n_star_from_angle(20.0);
    ASSERT_TRUE(nsr.message.empty());
    ASSERT_GE(nsr.n_star, 2);

    SplitParams params;
    params.n_star = nsr.n_star;
    auto result = split_segments(pslg, params);

    // 4 symmetric segments, each gets nsr.n_star subsegments
    EXPECT_EQ(result.all_segments.size(), 4u * nsr.n_star);
}

TEST(SegmentSplitter, AngleAbove30FallsBackToOne) {
    auto pslg = make_unit_square();

    auto nsr = compute_n_star_from_angle(33.0);
    EXPECT_EQ(nsr.n_star, 1);
    EXPECT_FALSE(nsr.message.empty());

    SplitParams params;
    params.n_star = nsr.n_star;
    auto result = split_segments(pslg, params);

    // n_star=1 means no extra splitting on symmetric geometry
    EXPECT_EQ(result.all_segments.size(), 4u);
}

// --- Geometry behavior tests ---

TEST(SegmentSplitter, ParallelFarApartFewerSplits) {
    PSLG far_pslg;
    far_pslg.vertices = {{0, 0}, {10, 0}, {0, 100}, {10, 100}};
    far_pslg.segments = {{0, 1}, {2, 3}};

    PSLG close_pslg;
    close_pslg.vertices = {{0, 0}, {10, 0}, {0, 0.1}, {10, 0.1}};
    close_pslg.segments = {{0, 1}, {2, 3}};

    SplitParams params;
    params.n_star = 3;

    auto far_result = split_segments(far_pslg, params);
    auto close_result = split_segments(close_pslg, params);

    EXPECT_LE(far_result.all_segments.size(), close_result.all_segments.size());
}

TEST(SegmentSplitter, LShapeMoreSplitsNearCorner) {
    PSLG pslg;
    pslg.vertices = {
        {0, 0}, {10, 0}, {10, 1}, {1, 1}, {1, 10}, {0, 10}
    };
    pslg.segments = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}
    };

    SplitParams params;
    params.n_star = 3;
    auto result = split_segments(pslg, params);

    EXPECT_GT(result.all_segments.size(), 6u);
}

TEST(SegmentSplitter, SingleSegmentNoSplitsDefaultNStar) {
    // With n_star=1, a single segment still gets just 1 subsegment.
    PSLG pslg;
    pslg.vertices = {{0, 0}, {5, 0}};
    pslg.segments = {{0, 1}};

    auto result = split_segments(pslg);
    EXPECT_EQ(result.all_segments.size(), 1u);
    EXPECT_EQ(result.all_vertices.size(), 2u);
}

TEST(SegmentSplitter, SplitDeterministic) {
    // Split the same PSLG twice, verify all vertex coordinates are bitwise identical
    PSLG pslg;
    pslg.vertices = {
        {0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}
    };
    pslg.segments = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}
    };

    SplitParams params;
    params.n_star = 5;
    auto r1 = split_segments(pslg, params);
    auto r2 = split_segments(pslg, params);

    ASSERT_EQ(r1.all_vertices.size(), r2.all_vertices.size());
    ASSERT_EQ(r1.all_segments.size(), r2.all_segments.size());

    for (size_t i = 0; i < r1.all_vertices.size(); ++i) {
        EXPECT_EQ(r1.all_vertices[i].x, r2.all_vertices[i].x)
            << "Vertex " << i << " x differs";
        EXPECT_EQ(r1.all_vertices[i].y, r2.all_vertices[i].y)
            << "Vertex " << i << " y differs";
    }
}

// --- Quadtree equivalence tests ---

namespace {

void assert_quadtree_equivalent(const PSLG& pslg, int n_star) {
    SplitParams naive_params;
    naive_params.n_star = n_star;
    naive_params.use_quadtree = false;
    auto naive = split_segments(pslg, naive_params);

    SplitParams qt_params;
    qt_params.n_star = n_star;
    qt_params.use_quadtree = true;
    auto qt = split_segments(pslg, qt_params);

    ASSERT_EQ(naive.all_vertices.size(), qt.all_vertices.size())
        << "Vertex count mismatch for n_star=" << n_star;
    ASSERT_EQ(naive.all_segments.size(), qt.all_segments.size())
        << "Segment count mismatch for n_star=" << n_star;

    for (size_t i = 0; i < naive.all_vertices.size(); ++i) {
        EXPECT_EQ(naive.all_vertices[i].x, qt.all_vertices[i].x)
            << "Vertex " << i << " x differs for n_star=" << n_star;
        EXPECT_EQ(naive.all_vertices[i].y, qt.all_vertices[i].y)
            << "Vertex " << i << " y differs for n_star=" << n_star;
    }

    for (size_t i = 0; i < naive.all_segments.size(); ++i) {
        EXPECT_EQ(naive.all_segments[i].p, qt.all_segments[i].p)
            << "Segment " << i << " p differs for n_star=" << n_star;
        EXPECT_EQ(naive.all_segments[i].q, qt.all_segments[i].q)
            << "Segment " << i << " q differs for n_star=" << n_star;
    }
}

} // namespace

TEST(QuadtreeEquivalence, UnitSquare) {
    PSLG pslg = make_unit_square();
    assert_quadtree_equivalent(pslg, 1);
    assert_quadtree_equivalent(pslg, 4);
    assert_quadtree_equivalent(pslg, 10);
}

TEST(QuadtreeEquivalence, LShape) {
    PSLG pslg;
    pslg.vertices = {
        {0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}
    };
    pslg.segments = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}
    };
    assert_quadtree_equivalent(pslg, 3);
    assert_quadtree_equivalent(pslg, 5);
}

TEST(QuadtreeEquivalence, NarrowChannel) {
    PSLG pslg;
    pslg.vertices = {
        {0, 0}, {10, 0}, {10, 0.2}, {0, 0.2}
    };
    pslg.segments = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}
    };
    assert_quadtree_equivalent(pslg, 3);
}

TEST(QuadtreeEquivalence, IsolatedSegment) {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {5, 0}};
    pslg.segments = {{0, 1}};
    assert_quadtree_equivalent(pslg, 1);
    assert_quadtree_equivalent(pslg, 10);
}

TEST(QuadtreeEquivalence, FarApartParallels) {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {10, 0}, {0, 100}, {10, 100}};
    pslg.segments = {{0, 1}, {2, 3}};
    assert_quadtree_equivalent(pslg, 3);
}

TEST(QuadtreeEquivalence, EmptyPSLG) {
    PSLG pslg;
    assert_quadtree_equivalent(pslg, 1);
}

TEST(QuadtreeEquivalence, StrayVertex) {
    PSLG pslg;
    pslg.vertices = {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}};
    pslg.segments = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    assert_quadtree_equivalent(pslg, 4);
}

TEST(QuadtreeEquivalence, LargeGrid) {
    // 20x20 grid with horizontal + vertical segments
    PSLG pslg;
    int grid = 20;
    // Create vertices: (grid+1) x (grid+1)
    for (int y = 0; y <= grid; ++y) {
        for (int x = 0; x <= grid; ++x) {
            pslg.vertices.push_back({static_cast<double>(x),
                                     static_cast<double>(y)});
        }
    }
    // Horizontal segments
    for (int y = 0; y <= grid; ++y) {
        for (int x = 0; x < grid; ++x) {
            int p = y * (grid + 1) + x;
            int q = p + 1;
            pslg.segments.push_back({p, q});
        }
    }
    // Vertical segments
    for (int x = 0; x <= grid; ++x) {
        for (int y = 0; y < grid; ++y) {
            int p = y * (grid + 1) + x;
            int q = p + (grid + 1);
            pslg.segments.push_back({p, q});
        }
    }
    assert_quadtree_equivalent(pslg, 3);
}

// --- Multi-threaded determinism tests ---

#ifdef _OPENMP
#include <omp.h>

namespace {

void assert_threaded_determinism(const PSLG& pslg, const SplitParams& base_params) {
    // Run with 1 thread as baseline
    SplitParams params1 = base_params;
    params1.num_threads = 1;
    auto baseline = split_segments(pslg, params1);

    for (int threads : {2, 4}) {
        SplitParams params_mt = base_params;
        params_mt.num_threads = threads;
        auto result = split_segments(pslg, params_mt);

        ASSERT_EQ(baseline.all_vertices.size(), result.all_vertices.size())
            << "Vertex count mismatch with " << threads << " threads";
        ASSERT_EQ(baseline.all_segments.size(), result.all_segments.size())
            << "Segment count mismatch with " << threads << " threads";

        for (size_t i = 0; i < baseline.all_vertices.size(); ++i) {
            EXPECT_EQ(baseline.all_vertices[i].x, result.all_vertices[i].x)
                << "Vertex " << i << " x differs with " << threads << " threads";
            EXPECT_EQ(baseline.all_vertices[i].y, result.all_vertices[i].y)
                << "Vertex " << i << " y differs with " << threads << " threads";
        }

        for (size_t i = 0; i < baseline.all_segments.size(); ++i) {
            EXPECT_EQ(baseline.all_segments[i].p, result.all_segments[i].p)
                << "Segment " << i << " p differs with " << threads << " threads";
            EXPECT_EQ(baseline.all_segments[i].q, result.all_segments[i].q)
                << "Segment " << i << " q differs with " << threads << " threads";
        }
    }
}

} // namespace

TEST(ThreadSafety, OpenMPRuntimeHonorsThreadCount) {
    int old_max = omp_get_max_threads();

    for (int requested : {2, 4}) {
        omp_set_num_threads(requested);
        int actual = 0;
        #pragma omp parallel
        {
            #pragma omp single
            {
                actual = omp_get_num_threads();
            }
        }
        EXPECT_EQ(actual, requested)
            << "OpenMP did not spawn " << requested << " threads; "
               "ThreadSafety tests cannot verify multi-threaded behavior";
    }

    omp_set_num_threads(old_max);
}

TEST(ThreadSafety, NaiveSquare) {
    PSLG pslg = make_unit_square();
    SplitParams params;
    params.n_star = 4;
    params.use_quadtree = false;
    assert_threaded_determinism(pslg, params);
}

TEST(ThreadSafety, NaiveLShape) {
    PSLG pslg;
    pslg.vertices = {
        {0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}
    };
    pslg.segments = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}
    };
    SplitParams params;
    params.n_star = 3;
    params.use_quadtree = false;
    assert_threaded_determinism(pslg, params);
}

TEST(ThreadSafety, QuadtreeSquare) {
    PSLG pslg = make_unit_square();
    SplitParams params;
    params.n_star = 4;
    params.use_quadtree = true;
    assert_threaded_determinism(pslg, params);
}

TEST(ThreadSafety, QuadtreeLShape) {
    PSLG pslg;
    pslg.vertices = {
        {0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}
    };
    pslg.segments = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}
    };
    SplitParams params;
    params.n_star = 3;
    params.use_quadtree = true;
    assert_threaded_determinism(pslg, params);
}

TEST(ThreadSafety, QuadtreeLargeGrid) {
    PSLG pslg;
    int grid = 20;
    for (int y = 0; y <= grid; ++y) {
        for (int x = 0; x <= grid; ++x) {
            pslg.vertices.push_back({static_cast<double>(x),
                                     static_cast<double>(y)});
        }
    }
    for (int y = 0; y <= grid; ++y) {
        for (int x = 0; x < grid; ++x) {
            int p = y * (grid + 1) + x;
            int q = p + 1;
            pslg.segments.push_back({p, q});
        }
    }
    for (int x = 0; x <= grid; ++x) {
        for (int y = 0; y < grid; ++y) {
            int p = y * (grid + 1) + x;
            int q = p + (grid + 1);
            pslg.segments.push_back({p, q});
        }
    }
    SplitParams params;
    params.n_star = 3;
    params.use_quadtree = true;
    assert_threaded_determinism(pslg, params);
}

#endif // _OPENMP

TEST(SegmentSplitter, FarthestEndpointCapsRefinement) {
    // Two parallel segments far apart. Without the farthest-endpoint cap,
    // F(x) would be large (distance to opposite segment) giving a large T.
    // With the cap, F(x) is clamped by max(x, l-x), limiting T.
    // The cap should reduce the total number of subsegments compared to
    // what we'd get if F were only the inter-segment distance.
    PSLG pslg;
    pslg.vertices = {{0, 0}, {10, 0}, {0, 100}, {10, 100}};
    pslg.segments = {{0, 1}, {2, 3}};

    SplitParams params;
    params.n_star = 3;
    auto result = split_segments(pslg, params);

    // With segments 100 apart and length 10, without the cap the inter-segment
    // distance (100) would dominate and T would be huge. With the cap, T is
    // bounded by the farthest-endpoint integral. Each segment should get
    // exactly n_star subsegments since both have the same T.
    EXPECT_EQ(result.all_segments.size(), 6u);  // 2 segments * 3 subsegments each
}
