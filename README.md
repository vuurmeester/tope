# tope - convex hull library in arbitrary dimension

A _polytope_ is the $d$ dimensional equivalent of a convex polygon.
Given $n$ points in $d$ dimensional Euclidian space (a _point cloud_), the convex hull is the smallest polytope containing all points.
Convex hulls can be used to create _Delaunay triangulations_, _Voronoi diagrams_ and _halfspace intersections_.

The standard 'go-to' for these types of problems is the qhull program <http://www.qhull.org>.
A small executable 'poly' is provided which works in a similar way to qhull.
In particular, 'poly' can parse 'rbox' output which facilitates comparison between the two programs.


# Some comparison

For example, creating a unit box, `rbox c | qhull` would give:

    Convex hull of 8 points in 3-d:

      Number of vertices: 8
      Number of facets: 6
      Number of non-simplicial facets: 6

    Statistics for: rbox c | qhull

      Number of points processed: 8
      Number of hyperplanes created: 11
      Number of distance tests for qhull: 34
      Number of distance tests for merging: 90
      Number of distance tests for checking: 48
      Number of merged facets: 6
      CPU seconds to compute hull (after input): 0.000154

While `rbox c | poly -m` gives:

    dimension 3
    npoints 8
      -0.5  -0.5  -0.5
      -0.5  -0.5   0.5
      -0.5   0.5  -0.5
      -0.5   0.5   0.5
       0.5  -0.5  -0.5
       0.5  -0.5   0.5
       0.5   0.5  -0.5
       0.5   0.5   0.5
    facets      = 6
    ridges      = 12
    verts       = 8
    time        = 3.91006e-05
    memory      = 4096

Note that you have to explicitly tell poly to attempt facet merging with the `-m` switch (this is default behavior in qhull).

For a more demanding calculation, we do 1000 points in 8 dimensions:
`rbox 1000 D8 | qhull` gives:

    Convex hull of 1000 points in 8-d:

      Number of vertices: 767
      Number of facets: 1309651
      Number of non-simplicial facets: 16

    Statistics for: rbox 1000 D8 | qhull

      Number of points processed: 846
      Number of hyperplanes created: 7968290
      Number of distance tests for qhull: 18205216
      Number of distance tests for merging: 71843039
      Number of distance tests for checking: 16617930
      Number of merged facets: 231
      CPU seconds to compute hull (after input): 39.22

Computing the same hull `rbox 1000 D8 | poly`:

    dimension 8
    npoints 1000
    facets      = 1309702
    ridges      = 5238808
    verts       = 767
    time        = 6.49203
    memory      = 1016463360

This happens to be one particular case where poly performs significantly faster than qhull.
Also, we are not trying to merge facets here.

