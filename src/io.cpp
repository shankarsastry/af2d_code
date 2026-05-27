#include "io.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace af2d {

namespace {

// Read a line, skipping blanks and comments (lines starting with #).
// Returns false on EOF.
bool read_line(std::ifstream& in, std::string& line) {
    while (std::getline(in, line)) {
        // Skip blank lines
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
        // Skip comment lines
        auto pos = line.find_first_not_of(" \t");
        if (pos != std::string::npos && line[pos] == '#') continue;
        return true;
    }
    return false;
}

// Strip inline comments: everything after # is removed.
std::string strip_comment(const std::string& line) {
    auto pos = line.find('#');
    if (pos != std::string::npos) return line.substr(0, pos);
    return line;
}

// Replace .poly or .node extension with .node
std::string companion_node_path(const std::string& poly_path) {
    auto dot = poly_path.rfind('.');
    if (dot == std::string::npos) return poly_path + ".node";
    return poly_path.substr(0, dot) + ".node";
}

// Parse vertices from lines into pslg. Returns the index offset (0 or 1).
int parse_vertices(std::ifstream& in, PSLG& pslg, int num_vertices,
                   int num_attributes, int has_boundary_markers) {
    int index_offset = 0;
    for (int i = 0; i < num_vertices; ++i) {
        std::string line;
        if (!read_line(in, line)) {
            throw std::runtime_error("Unexpected EOF reading vertices");
        }
        std::istringstream iss(strip_comment(line));
        int id;
        double x, y;
        iss >> id >> x >> y;
        if (iss.fail()) {
            throw std::runtime_error("Failed to parse vertex line: " + line);
        }
        if (i == 0) {
            index_offset = id;  // detect 0-based vs 1-based
        }
        // Skip attributes and boundary marker
        pslg.vertices.push_back({x, y});
    }
    return index_offset;
}

} // anonymous namespace

PSLG read_node(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open file: " + path);
    }

    std::string line;
    if (!read_line(in, line)) {
        throw std::runtime_error("Empty .node file: " + path);
    }

    std::istringstream header(strip_comment(line));
    int num_vertices, dim, num_attributes, has_boundary_markers;
    header >> num_vertices >> dim >> num_attributes >> has_boundary_markers;
    if (header.fail() || dim != 2) {
        throw std::runtime_error("Invalid .node header: " + line);
    }

    PSLG pslg;
    parse_vertices(in, pslg, num_vertices, num_attributes, has_boundary_markers);
    return pslg;
}

PSLG read_poly(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open file: " + path);
    }

    std::string line;

    // --- Vertex section ---
    if (!read_line(in, line)) {
        throw std::runtime_error("Empty .poly file: " + path);
    }

    std::istringstream vheader(strip_comment(line));
    int num_vertices, dim, num_attributes, has_boundary_markers;
    vheader >> num_vertices >> dim >> num_attributes >> has_boundary_markers;
    if (vheader.fail()) {
        throw std::runtime_error("Invalid .poly vertex header: " + line);
    }

    PSLG pslg;
    int index_offset = 0;

    if (num_vertices == 0) {
        // Read vertices from companion .node file.
        // Detect index base from the first segment line instead.
        std::string node_path = companion_node_path(path);
        pslg = read_node(node_path);
        index_offset = -1;  // sentinel: detect from first segment ID
    } else {
        if (dim != 2) {
            throw std::runtime_error("Invalid dimension in .poly header: " + line);
        }
        index_offset = parse_vertices(in, pslg, num_vertices, num_attributes,
                                      has_boundary_markers);
    }

    // --- Segment section ---
    if (!read_line(in, line)) {
        return pslg;  // No segments section
    }

    std::istringstream sheader(strip_comment(line));
    int num_segments, seg_boundary_markers;
    sheader >> num_segments >> seg_boundary_markers;
    if (sheader.fail()) {
        throw std::runtime_error("Invalid .poly segment header: " + line);
    }

    for (int i = 0; i < num_segments; ++i) {
        if (!read_line(in, line)) {
            throw std::runtime_error("Unexpected EOF reading segments");
        }
        std::istringstream iss(strip_comment(line));
        int id, p, q;
        iss >> id >> p >> q;
        if (iss.fail()) {
            throw std::runtime_error("Failed to parse segment line: " + line);
        }
        // When vertices came from a companion .node (index_offset == -1),
        // detect the index base from the smallest vertex reference we see.
        if (index_offset == -1) {
            index_offset = std::min(p, q);  // 0 or 1
        }
        pslg.segments.push_back({p - index_offset, q - index_offset});
    }

    // --- Holes section ---
    if (!read_line(in, line)) {
        return pslg;  // No holes section
    }

    std::istringstream hheader(strip_comment(line));
    int num_holes;
    hheader >> num_holes;
    if (hheader.fail()) {
        throw std::runtime_error("Invalid .poly holes header: " + line);
    }

    for (int i = 0; i < num_holes; ++i) {
        if (!read_line(in, line)) {
            throw std::runtime_error("Unexpected EOF reading holes");
        }
        std::istringstream iss(strip_comment(line));
        int id;
        double x, y;
        iss >> id >> x >> y;
        if (iss.fail()) {
            throw std::runtime_error("Failed to parse hole line: " + line);
        }
        pslg.holes.push_back({x, y});
    }

    return pslg;
}

void write_node(const std::string& path, const TriangleOutput& out) {
    std::ofstream f(path);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + path);
    }

    f << out.num_points << " 2 0 0\n";
    f.precision(17);
    for (int i = 0; i < out.num_points; ++i) {
        f << i << " " << out.points[i].x << " " << out.points[i].y << "\n";
    }
}

void write_ele(const std::string& path, const TriangleOutput& out) {
    std::ofstream f(path);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + path);
    }

    f << out.num_triangles << " 3 0\n";
    for (int i = 0; i < out.num_triangles; ++i) {
        f << i << " " << out.triangles[i][0]
          << " " << out.triangles[i][1]
          << " " << out.triangles[i][2] << "\n";
    }
}

void write_edge(const std::string& path, const TriangleOutput& out) {
    std::ofstream f(path);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + path);
    }

    f << out.num_edges << " 0\n";
    for (int i = 0; i < out.num_edges; ++i) {
        f << i << " " << out.edges[i].p << " " << out.edges[i].q << "\n";
    }
}

} // namespace af2d
