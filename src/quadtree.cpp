#include "quadtree.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace af2d {

// --- BBox ---

bool BBox::contains(double x, double y) const {
    return x >= x_min && x <= x_max && y >= y_min && y <= y_max;
}

bool BBox::intersects(const BBox& other) const {
    return x_min <= other.x_max && x_max >= other.x_min &&
           y_min <= other.y_max && y_max >= other.y_min;
}

BBox BBox::expanded(double pad) const {
    return {x_min - pad, y_min - pad, x_max + pad, y_max + pad};
}

BBox BBox::from_segment(const Point2D& a, const Point2D& b) {
    return {std::min(a.x, b.x), std::min(a.y, b.y),
            std::max(a.x, b.x), std::max(a.y, b.y)};
}

// --- Quadtree::Node ---

bool Quadtree::Node::is_leaf() const {
    return !children[0];
}

// --- Quadtree ---

void Quadtree::build(const PSLG& pslg) {
    if (pslg.vertices.empty()) {
        root_.reset();
        return;
    }

    // Compute bounding box of all vertices
    double xmin = pslg.vertices[0].x, xmax = xmin;
    double ymin = pslg.vertices[0].y, ymax = ymin;
    for (const auto& v : pslg.vertices) {
        xmin = std::min(xmin, v.x);
        xmax = std::max(xmax, v.x);
        ymin = std::min(ymin, v.y);
        ymax = std::max(ymax, v.y);
    }

    // Add small padding to avoid degenerate zero-size bounding box
    double pad = std::max(1e-10, std::max(xmax - xmin, ymax - ymin) * 1e-6);
    root_ = std::make_unique<Node>();
    root_->bounds = {xmin - pad, ymin - pad, xmax + pad, ymax + pad};

    // Insert all vertices
    for (int i = 0; i < static_cast<int>(pslg.vertices.size()); ++i) {
        const auto& v = pslg.vertices[i];
        FeatureItem item{FeatureType::Vertex, i, {v.x, v.y, v.x, v.y}};
        insert(root_.get(), item, 0);
    }

    // Insert all segments
    for (int i = 0; i < static_cast<int>(pslg.segments.size()); ++i) {
        const auto& s = pslg.segments[i];
        BBox bbox = BBox::from_segment(pslg.vertices[s.p], pslg.vertices[s.q]);
        FeatureItem item{FeatureType::Segment, i, bbox};
        insert(root_.get(), item, 0);
    }
}

void Quadtree::insert(Node* node, const FeatureItem& item, int depth) {
    if (!node->bounds.intersects(item.bbox)) {
        return;
    }

    if (node->is_leaf()) {
        node->items.push_back(item);
        if (static_cast<int>(node->items.size()) > MAX_ITEMS_PER_NODE &&
            depth < MAX_DEPTH) {
            subdivide(node);
        }
        return;
    }

    // Try to insert into children
    bool inserted = false;
    for (auto& child : node->children) {
        if (child->bounds.intersects(item.bbox)) {
            // Check if item fits entirely within this child
            if (item.bbox.x_min >= child->bounds.x_min &&
                item.bbox.x_max <= child->bounds.x_max &&
                item.bbox.y_min >= child->bounds.y_min &&
                item.bbox.y_max <= child->bounds.y_max) {
                insert(child.get(), item, depth + 1);
                return;
            }
            inserted = true;
        }
    }

    // Item spans multiple children — keep in parent
    if (inserted) {
        node->items.push_back(item);
    }
}

void Quadtree::subdivide(Node* node) {
    double mx = (node->bounds.x_min + node->bounds.x_max) * 0.5;
    double my = (node->bounds.y_min + node->bounds.y_max) * 0.5;

    node->children[0] = std::make_unique<Node>();
    node->children[0]->bounds = {node->bounds.x_min, node->bounds.y_min, mx, my};

    node->children[1] = std::make_unique<Node>();
    node->children[1]->bounds = {mx, node->bounds.y_min, node->bounds.x_max, my};

    node->children[2] = std::make_unique<Node>();
    node->children[2]->bounds = {node->bounds.x_min, my, mx, node->bounds.y_max};

    node->children[3] = std::make_unique<Node>();
    node->children[3]->bounds = {mx, my, node->bounds.x_max, node->bounds.y_max};

    // Redistribute items
    std::vector<FeatureItem> old_items = std::move(node->items);
    node->items.clear();

    for (const auto& item : old_items) {
        // Check if item fits entirely in one child
        bool placed = false;
        for (auto& child : node->children) {
            if (item.bbox.x_min >= child->bounds.x_min &&
                item.bbox.x_max <= child->bounds.x_max &&
                item.bbox.y_min >= child->bounds.y_min &&
                item.bbox.y_max <= child->bounds.y_max) {
                child->items.push_back(item);
                placed = true;
                break;
            }
        }
        if (!placed) {
            // Spans multiple children — stays in parent
            node->items.push_back(item);
        }
    }
}

void Quadtree::query(const BBox& region, std::vector<FeatureItem>& results) const {
    if (!root_) return;
    query_impl(root_.get(), region, results);
}

void Quadtree::query_impl(const Node* node, const BBox& region,
                           std::vector<FeatureItem>& results) const {
    if (!node->bounds.intersects(region)) {
        return;
    }

    for (const auto& item : node->items) {
        if (item.bbox.intersects(region)) {
            results.push_back(item);
        }
    }

    if (!node->is_leaf()) {
        for (const auto& child : node->children) {
            query_impl(child.get(), region, results);
        }
    }
}

// --- Distance helpers ---

double min_distance_point_to_segment(const Point2D& p,
                                     const Point2D& a, const Point2D& b) {
    Point2D ab = b - a;
    double ab_len_sq = ab.dot(ab);
    if (ab_len_sq < 1e-30) {
        return p.distance(a);
    }
    double t = (p - a).dot(ab) / ab_len_sq;
    t = std::clamp(t, 0.0, 1.0);
    Point2D proj = a + ab * t;
    return p.distance(proj);
}

double min_distance_segment_to_segment(const Point2D& a1, const Point2D& b1,
                                       const Point2D& a2, const Point2D& b2) {
    // Minimum of the 4 endpoint-to-segment distances
    double d1 = min_distance_point_to_segment(a1, a2, b2);
    double d2 = min_distance_point_to_segment(b1, a2, b2);
    double d3 = min_distance_point_to_segment(a2, a1, b1);
    double d4 = min_distance_point_to_segment(b2, a1, b1);
    return std::min({d1, d2, d3, d4});
}

} // namespace af2d
