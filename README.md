# AF2D

Adaptive 2D mesh generator using LFS-proportional segment splitting and constrained Delaunay triangulation. Implements the algorithm from [Sastry (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8078854/).

Accepts Triangle's `.node` and `.poly` file formats as input and produces `.node`, `.ele`, and `.edge` files as output, making it a drop-in replacement for Triangle's workflow.

## Building

Requires CMake 3.24+ and a C++20 compiler. Dependencies (Triangle, GoogleTest) are fetched automatically.

```
cmake -B build
cmake --build build
```

## Running

```
af2d_main <input.poly|input.node> [-q<angle>] [-n<n_star>] [-t]
```

### Arguments

**Positional:**

- `<input>` -- Input geometry file. `.poly` files provide vertices, segments, and holes. `.node` files provide vertices only (Triangle will compute the convex hull).

**Options:**

- `-q<angle>` -- Minimum angle in degrees. This does two things:
  1. If the angle is below 30, computes n\* from the CDT non-encroachment condition (Lemma 4.12, Sastry 2021) to determine segment splitting density.
  2. Passes the angle constraint to Triangle, which inserts Steiner points until no triangle has an angle below the threshold.

  If the angle is 30 or above, n\* defaults to 1 (the formula has a singularity), but Triangle still enforces the angle constraint via refinement.

- `-n<n_star>` -- Explicit n\* value for segment splitting. Controls how many subsegments the shortest reference segment gets; other segments get proportionally more based on their local feature size ratio. Overrides the n\* derived from `-q` if both are given. Does not affect Triangle's angle refinement.

- `-t` -- Enable quadtree acceleration for LFS computation. Builds a spatial index to prune distant features during segment splitting, reducing the cost of the pairwise distance computations. Produces bitwise-identical output to the default algorithm. This feature has been unit tested for correctness (including equivalence tests against the naive path on several geometries), but has not been stress tested on large real-world inputs.

### Precedence

When both `-q` and `-n` are given, `-n` sets n\* for segment splitting, and `-q` is still passed to Triangle for angle enforcement.

### Examples

```
# Default: n*=1, no angle refinement
af2d_main mesh.poly

# Minimum angle of 20 degrees (computes n*=13 from theory, Triangle enforces 20)
af2d_main mesh.poly -q20

# Minimum angle of 32 degrees (n* defaults to 1, Triangle enforces 32)
af2d_main mesh.poly -q32

# Explicit n*=4, no angle refinement
af2d_main mesh.poly -n4

# Explicit n*=4 with 25-degree angle enforcement
af2d_main mesh.poly -n4 -q25

# Quadtree-accelerated splitting
af2d_main mesh.poly -n4 -t
```

### Output

Output files are written alongside the input with `.1` appended to the basename:

- `<basename>.1.node` -- Vertex coordinates
- `<basename>.1.ele` -- Triangle connectivity (3 vertices per triangle)
- `<basename>.1.edge` -- Edge list

All output uses 0-based indexing.

A summary is printed to stdout:

```
Input:     6 vertices, 6 segments, 0 holes
n*:        4
Output:    42 vertices, 56 triangles, 97 edges
Files:     mesh.1.node
           mesh.1.ele
           mesh.1.edge
```

## Input file formats

AF2D uses Triangle's standard file formats.

### .node (vertices)

```
<num_vertices> 2 <num_attributes> <0|1 boundary markers>
<id> <x> <y> [attributes...] [boundary_marker]
...
```

### .poly (vertices + segments + holes)

```
<num_vertices> 2 <num_attributes> <0|1 boundary markers>
<id> <x> <y> [attributes...] [boundary_marker]
...
<num_segments> <0|1 boundary markers>
<id> <endpoint1> <endpoint2> [boundary_marker]
...
<num_holes>
<id> <x> <y>
...
```

If the vertex count in the `.poly` header is 0, vertices are read from a companion `.node` file with the same basename. Both 0-based and 1-based indexing are detected automatically. Comments (`#`) and blank lines are skipped.

## Running tests

```
ctest --test-dir build --output-on-failure
```

## Algorithm overview

1. **Read** input PSLG (vertices, segments, holes)
2. **Compute local feature size (LFS)** for each segment as the lower envelope of distances to non-adjacent vertices, non-adjacent segments, and the segment's own farthest endpoint
3. **Solve the ODE** M'(t) = F(M(t)) where F is the LFS function, mapping each physical segment to a reference segment of length T
4. **Split segments** uniformly in reference space: n_i = floor(n\* * T_i / T_min) subsegments per segment, placed at positions M(t_k) in physical space
5. **Triangulate** the split PSLG using Shewchuk's Triangle library with constrained Delaunay triangulation
6. **Write** output mesh files
