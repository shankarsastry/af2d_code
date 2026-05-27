#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

namespace af2d {

struct Point2D {
    double x = 0.0;
    double y = 0.0;

    Point2D operator+(const Point2D& o) const { return {x + o.x, y + o.y}; }
    Point2D operator-(const Point2D& o) const { return {x - o.x, y - o.y}; }
    Point2D operator*(double s) const { return {x * s, y * s}; }
    Point2D operator/(double s) const { return {x / s, y / s}; }

    double dot(const Point2D& o) const { return x * o.x + y * o.y; }
    double norm() const { return std::sqrt(x * x + y * y); }
    double distance(const Point2D& o) const { return (*this - o).norm(); }
};

inline Point2D operator*(double s, const Point2D& p) { return p * s; }

struct Segment {
    int p = -1;
    int q = -1;
};

class PSLG {
public:
    std::vector<Point2D> vertices;
    std::vector<Segment> segments;
    std::vector<Point2D> holes;

    double segment_length(int i) const;
    Point2D point_on_segment(int i, double x) const;
    bool is_adjacent(int seg_i, int feature_seg) const;
    bool is_adjacent_vertex(int seg_i, int vert_idx) const;
    void validate() const;
};

} // namespace af2d
