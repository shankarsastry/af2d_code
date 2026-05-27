#pragma once

#include "pslg.h"
#include "triangle_wrapper.h"
#include <string>

namespace af2d {

// Read a .node file -> PSLG with vertices only (no segments/holes).
PSLG read_node(const std::string& path);

// Read a .poly file -> full PSLG (vertices, segments, holes).
// If vertex count is 0 in the .poly header, reads vertices from companion .node file.
PSLG read_poly(const std::string& path);

// Write output mesh files in Triangle's format (0-based indexing).
void write_node(const std::string& path, const TriangleOutput& out);
void write_ele(const std::string& path, const TriangleOutput& out);
void write_edge(const std::string& path, const TriangleOutput& out);

} // namespace af2d
