#pragma once

#include "pslg.h"
#include <array>
#include <memory>
#include <vector>

namespace af2d {

struct BBox {
    double x_min, y_min, x_max, y_max;

    bool contains(double x, double y) const;
    bool intersects(const BBox& other) const;
    BBox expanded(double pad) const;
    static BBox from_segment(const Point2D& a, const Point2D& b);
};

enum class FeatureType { Vertex, Segment };

struct FeatureItem {
    FeatureType type;
    int index;    // vertex or segment index into PSLG
    BBox bbox;
};

class Quadtree {
public:
    void build(const PSLG& pslg);
    void query(const BBox& region, std::vector<FeatureItem>& results) const;

private:
    static constexpr int MAX_ITEMS_PER_NODE = 16;
    static constexpr int MAX_DEPTH = 20;

    struct Node {
        BBox bounds;
        std::vector<FeatureItem> items;
        std::array<std::unique_ptr<Node>, 4> children;
        bool is_leaf() const;
    };

    std::unique_ptr<Node> root_;

    void insert(Node* node, const FeatureItem& item, int depth);
    void subdivide(Node* node);
    void query_impl(const Node* node, const BBox& region,
                    std::vector<FeatureItem>& results) const;
};

// Distance helpers for quadtree pruning
double min_distance_point_to_segment(const Point2D& p,
                                     const Point2D& a, const Point2D& b);
double min_distance_segment_to_segment(const Point2D& a1, const Point2D& b1,
                                       const Point2D& a2, const Point2D& b2);

} // namespace af2d
