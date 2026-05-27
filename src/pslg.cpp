#include "pslg.h"

namespace af2d {

double PSLG::segment_length(int i) const {
    const auto& s = segments.at(i);
    return vertices.at(s.p).distance(vertices.at(s.q));
}

Point2D PSLG::point_on_segment(int i, double x) const {
    const auto& s = segments.at(i);
    const auto& pa = vertices.at(s.p);
    const auto& pb = vertices.at(s.q);
    double l = pa.distance(pb);
    if (l < 1e-15) return pa;
    double t = x / l;
    return pa + (pb - pa) * t;
}

bool PSLG::is_adjacent(int seg_i, int feature_seg) const {
    const auto& si = segments.at(seg_i);
    const auto& sf = segments.at(feature_seg);
    return si.p == sf.p || si.p == sf.q || si.q == sf.p || si.q == sf.q;
}

bool PSLG::is_adjacent_vertex(int seg_i, int vert_idx) const {
    const auto& si = segments.at(seg_i);
    return si.p == vert_idx || si.q == vert_idx;
}

void PSLG::validate() const {
    int n = static_cast<int>(vertices.size());
    for (int i = 0; i < static_cast<int>(segments.size()); ++i) {
        const auto& s = segments[i];
        if (s.p < 0 || s.p >= n || s.q < 0 || s.q >= n) {
            throw std::invalid_argument("Segment index out of range");
        }
        if (s.p == s.q) {
            throw std::invalid_argument("Degenerate segment (same endpoints)");
        }
        double l = vertices[s.p].distance(vertices[s.q]);
        if (l < 1e-15) {
            throw std::invalid_argument("Zero-length segment");
        }
    }
}

} // namespace af2d
