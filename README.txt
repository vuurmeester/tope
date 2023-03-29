polytoop - convex hull library in arbitrary dimension

Given n points in d dimensional Euclidian space, the convex hull is the
smallest polytope containing all points. Convex hulls are used to create
Delaunay triangulations, Voronoi diagrams and halfspace intersections.

This library provides an alternative to qhull <http://www.qhull.org>, which
is a great resource for learning about the underlying math and implementation
details.
  - readability / extensibility
  - library interface
  - handling of numerical roundoff (e.g. infinitesimal facets, 'joggling', etc)
