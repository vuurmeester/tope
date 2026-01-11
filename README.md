# tope - convex hull library in arbitrary dimension

Given $n$ points in $d$ dimensional Euclidian space, the convex hull is the smallest polytope containing all points.
Convex hulls can be used to create Delaunay triangulations, Voronoi diagrams and halfspace intersections.

The standard 'go-to' for these types of problems is the qhull program <http://www.qhull.org>.
A small executable 'poly' is provided which works in a similar way to qhull.
In particular, 'poly' can parse 'rbox' output which facilitates comparison between the two programs.
