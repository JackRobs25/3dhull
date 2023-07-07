# hull3d


Our project implements the incremental 3d hull algorithm. This starts by creating a hull by making four faces using the first 4 points. Then it iterates through
the remaining points. For each point it determines if the point is inside the current hull, if it is it can be ignored, if it is outside we remove all the faces of the current hull that are visible from the point and reconnect the hull to this new point. This process repeats until all points have been considered, at which point the hull is complete and all faces are extreme.

In order to run the code you simply provide hull3d with the number of points you wish to plot. Once you have run the hull3d program you can press 'c' to fill the hull faces with color and 'a' to animate the hull process so you can see the hull being built "incrementally".

list of main functions:

// returns true if point p is inside the shape defined by the faces in hull
bool point_inside_convex_polygon(point3d p, vector<triangle3d>& hull);

// return true if edge is shared with any of the edges of face
bool shares_edge(pair<point3d, point3d> edge, triangle3d face)

// return true if edge1 < edge 2. Ordered by x, then y, then z
bool edgesCmp(edge3d edge1, edge3d edge2)

//orient the first face containing points 0,1,2 so that it faces away from point 3
void orient_first_face(triangle3d* t, point3d p)

// orient so that the normal is outside the hull 
// ie all the points in the hull are behind the face or coplanar with it
void orient_triangle(triangle3d* t, vector<triangle3d>& hull)

/* compute the convex hull of the points */
void incremental_hull(vector<point3d>& points, vector<triangle3d>& hull)