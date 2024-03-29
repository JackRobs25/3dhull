#ifndef __geom_h
#define __geom_h


#include <vector>

using namespace std; 




typedef struct _point2d {
  int  x,y; 
} point2d;



typedef struct _point3d {
  int  x,y,z; 
} point3d;



typedef struct _triangle3d {
  int ia, ib, ic; //indices in the array of points, technically not
		  //necessary, but makes writing the code easier and
		  //debugging easier
  point3d *a,*b,*c;   //to avoid duplication of data a triangle stores
		      //pointers to the points
  double color[3];  //this is the color used to render this triangle
} triangle3d;


typedef struct _edge3d {
  int ia, ib; //indices in the array of points, technically not
	      //necessary, but makes writing the code easier and
	      //debugging easier
  point3d *a, *b; //pointers to the points
} edge3d; 





/*  ****************** 2D geometric functions ****************** */

/* **************************************** */
/* returns 2 x the signed area of triangle abc */
long long  signed_area2d(point2d a, point2d b, point2d c);

/* return true if p,q,r collinear, and false  otherwise */
bool collinear(point2d p, point2d q, point2d r); 

/* return true if c is  strictly left of ab; false otherwise */
bool left_strictly(point2d a, point2d b, point2d c); 

/* return True if c is   left of ab or on ab; false otherwise */
bool left_on(point2d a, point2d b, point2d c); 

/* return true if c is  strictly right  of ab; 0 otherwise */
bool right_strictly(point2d a, point2d b, point2d c); 

long long  dist2d(point2d a, point2d b);





/*  ****************** 3D geometric functions ****************** */

/* returns 6 times the signed volume of abcd. The volume is positive
   if d is behind abc (i.e. on opposite side as the normal); negative
   if d is in front (i.e. same side as the normal) of abc, and 0 if
   abcd are coplanar.
 */
long long signed_volume(point3d a, point3d b, point3d c, point3d d);


/* return True if points are on the same plane, and False otherwise */
bool coplanar(point3d a, point3d b, point3d c, point3d d);


/* return True if d is  strictly in front of abc; False otherwise */
bool infront (point3d a, point3d b, point3d c, point3d d); 


/* return true if face  defined by points (i,j,k) is extreme */
bool face_is_extreme(int i, int j, int k,  vector<point3d>& points);

void print_triangle(triangle3d t, const char* l1, const char* l2);

void print_faces(vector<triangle3d> hull);



/* compute the convex hull of the points */
void naive_hull(vector<point3d>& points, vector<triangle3d>& hull);
void incremental_hull(vector<point3d>& points, vector<triangle3d> & hull);
void gift_wrapping_hull(vector<point3d>& points, vector<triangle3d>& hull);


#endif
