#include "geom.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>

using namespace std; 






/* ************************************************************ */
/*  ****************** 2D functions ****************** */
/* ************************************************************ */


/* **************************************** */
/* returns 2 x the signed area of triangle abc */
long long signed_area2d(point2d a, point2d b, point2d c) {
  return  (long long)(2* ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)));
}

/* **************************************** */
/* return true  if a, b, c collinear, and false otherwise */
bool collinear(point2d a, point2d b, point2d c) {
  long long area = signed_area2d(a, b, c);
  //return  (area < EPSILON &&  area > -EPSILON); 
  return (area == 0);
}

/* **************************************** */
/* return True if c is  strictly left of ab; false otherwise */
bool  left_strictly(point2d a, point2d b, point2d c) {
  //return (signed_area2d(a, b, c) > EPSILON);
  return (signed_area2d(a, b, c) > 0);
}

/* return True if c is   left of ab or on ab; false otherwise */
bool  left_on(point2d a, point2d b, point2d c) {
  //return (signed_area2d(a, b, c) > -EPSILON);
  return (signed_area2d(a, b, c) >= 0);
}


/* **************************************** */
/* return true if c is  strictly right  of ab; false otherwise */
bool  right_strictly(point2d a, point2d b, point2d c) {
  //return (signed_area2d(a, b, c) < -EPSILON);
  return (signed_area2d(a, b, c) < 0);
}

/* **************************************** */
long long  dist2d(point2d a, point2d b) {
  return (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y);
}







//////////////////////////////////////////////////////////////////////
/*  ****************** 3D functions ****************** */
//////////////////////////////////////////////////////////////////////



long long  det2(long   a1, long   a2,  long  b1,  long  b2) {
  return a1*b2 - a2*b1; 
}

long long  det3( long  a1,  long  a2,  long  a3,
		 long  b1,  long  b2,  long  b3,
		 long   c1,  long  c2,  long  c3) {
  
  long long res;
  res =  (long long)a1* det2(b2, b3, c2, c3);
  res -= (long long) a2*det2(b1, b3, c1, c3);
  res += (long long) a3*det2(b1, b2, c1, c2);
  return res; 
}

/* ************************************************************ */
/* returns 6 times the signed volume of abcd. The volume is positive
   if d is behind abc (i.e. on opposite side as the normal); negative
   if d is in front (i.e. same side as the normal) of abc, and 0 if
   abcd are coplanar.
 */
long long  signed_volume(point3d a, point3d b, point3d c, point3d d) {
  long long res;
  res =  (long long)a.x*det3(b.y, b.z, 1, c.y, c.z, 1, d.y, d.z, 1); 
  res -= (long long) a.y*det3(b.x, b.z, 1, c.x, c.z, 1, d.x, d.z, 1);
  res += (long long)a.z*det3(b.x, b.y, 1, c.x, c.y, 1, d.x, d.y, 1); 
  res -= (long long)det3(b.x, b.y, b.z, c.x, c.y, c.z, d.x, d.y, d.z);
  return res; 
}

/* ************************************************************ */
/* return True if points are on the same plane, and False otherwise */
bool  coplanar(point3d a, point3d b, point3d c, point3d d) {

  long long  vol = signed_volume(a, b, c, d); 
  // return  (vol < EPSILON && vol > -EPSILON); //if using doubles
  return  (vol == 0); 
}

/* ************************************************************ */
/* return True if d is  strictly in front of abc; False otherwise */
bool  infront_strictly (point3d a, point3d b, point3d c, point3d d) {

  long long  vol = signed_volume(a, b, c, d); 
  // return  (vol < -EPSILON); //if using doubles
  return  (vol < 0); 
}

/* ************************************************************ */
//return true if face  defined by points (i,j,k) is extreme 
bool face_is_extreme(int i, int j, int k,  vector<point3d>& points) {

  int nfront=0, nback=0; 
  for (int l =0; l< points.size(); l++) {
    if (l==i || l==j || l==k)  continue;
    if (coplanar(points[i], points[j], points[k], points[l])) continue; 
    if (infront_strictly(points[i], points[j], points[k], points[l])) nfront++; 
    else  nback++; 
    
    if (nfront == 0 || nback == 0) return false; 
  }
  //if we got here, all on same side 
  return true;
}





/* compute the convex hull of the points */
void naive_hull(vector<point3d>& points, vector<triangle3d>& hull) {
  for (int i = 0; i < points.size(); i++){
    for (int j = i+1; j < points.size(); j++){
      for (int k = j+1; k < points.size(); k++){
        // check each unique face made up of three points 
        if (face_is_extreme(i, j, k, points)){
          triangle3d currFace; 
          currFace.a = &points[i];
          currFace.b = &points[j];
          currFace.c = &points[k];
          currFace.color[0] = (double)rand() / RAND_MAX; // R component
          currFace.color[1] = (double)rand() / RAND_MAX; // G component
          currFace.color[2] = (double)rand() / RAND_MAX; // B component
          hull.push_back(currFace);
        }
      }// for k
    }// for j
  }// for i
  /*for (int i = 0; i<hull.size(); i++){
    if(!face_is_extreme(hull[i].ia, hull[i].ib, hull[i].ic, points)){
      print_triangle(hull[i], "face", "not extreme\n");
      exit(1);
    }
  }*/
}

void print_triangle(triangle3d t, const char* l1, const char* l2){
  printf("%s", l1);
  printf("(%d, %d, %d)", t.ia, t.ib, t.ic);
  printf("(%d, %d, %d), (%d, %d, %d), (%d, %d, %d)", t.a->x, t.a->y, t.a->z, t.b->x, t.b->y, t.b->z, t.c->x, t.c->y, t.c->z);
  printf("%s", l2);
}

void print_faces(vector<triangle3d> hull){
  printf("hull faces: ");
  for (int i = 0; i<hull.size(); i++){
    print_triangle(hull[i], " ", " ");
  }
  printf("\n");
}



/*
//////////////////////////////////////////////////////////////////////

                               GIFT WRAPPING
//////////////////////////////////////////////////////////////////////
*/

//helper functions
void print_point(point3d p, int i) {
  printf("%3d(%d, %d, %d) ", i, p.x, p.y, p.z); 
}

/* ************************************************************ */
//return index of  point with max x-coord 
int find_right_most_point(vector<point3d>& points) {

  if(points.size() == 0) return -1;
  
  int rightmost = 0;
  for (int i=1; i< points.size(); i++){
    if (points[rightmost].x < points[i].x) {
      rightmost = i;
    }
  }
  return rightmost; 
}


  
/* ************************************************************ */
//return true if the edge between points i, j is extreme in the 2d
//projection of all points on the z-plane
bool is_edge_projection_extreme(int first_point, int second_point, vector<point3d>& points) {
  
  point2d first = {points[first_point].x, points[first_point].y};
  point2d second = {points[second_point].x, points[second_point].y}; 
  
  for (int i=0; i<points.size(); i++) {

    if (i==first_point || i==second_point) continue;
    
    point2d p = {points[i].x, points[i].y};
    if (collinear(first, second, p)) continue; 
    
    if (right_strictly(first, second, p)) {
      printf("ERROR first EDGE  NOT extreme\n"); 
      return false; 
    }
  }//for
  return true; 
}



/* ************************************************************ */
//return an edge on the hull 
edge3d find_first_edge_on_hull(vector<point3d>& points) {

  int first_index = find_right_most_point(points);
  printf("%15s", "first point: ");
  print_point(points[first_index], first_index);
  printf("\n"); 

  //project all points onto z=0 plane and 2d gift-wrap to find first edge from first_point
  point2d first = {points[first_index].x, points[first_index].y};
  int second_index = -1;
  point2d second = {0, 0}, p;
  for (int i=0; i<points.size(); i++) {

    if (i==first_index) continue;

    //current point 
    p.x = points[i].x; p.y= points[i].y;

    if (second_index==-1) {
      second_index = i;
      second.x =  points[i].x;
      second.y = points[i].y;
    } else {
      if (right_strictly(first, second, p)) {
	second_index = i;
	second.x = points[i].x;
	second.y = points[i].y;
      }
    }
  }//for
  
  printf("%15s", "second point: ");
  print_point(points[second_index], second_index);
  printf("\n"); 

  //sanity check that edge is indeed extreme 
  assert(is_edge_projection_extreme(first_index, second_index, points)); 
  
  edge3d e = {first_index, second_index, &points[first_index], &points[second_index] }; 
  return e; 
}




/* ************************************************************ */
//p, q are indices of two points, edge (p,q) assumed to be extreme ie on the hull
//returns the index r of  vertex which  is front-most as seen from pq 
int pivot_around_edge(int p, int  q, vector<point3d>& points) {

  printf("pivot around edge (%d, %d): \n", p, q); 

  int r = -1; 
  //find the front-most point
  for (int i=0; i< points.size(); i++) {
    if (i==p || i==q ) continue;
    
    if (r == -1) {
     r = i; 
    } else {
      if (infront_strictly(points[p],  points[q], points[r],  points[i])) {
	//this point is more infront than r
	//printf("\tr=%d. point %d in front  of %d, swapping\n", r, i, r); 
	r = i;
      }
    }//else 
  }//for 
  
  //sanity check
  if (!face_is_extreme(p, q, r, points)) {
    printf("pivot_around_edge: returns a face [%d,%d,%d] that's NOT EXTREME\n", p, q, r);
    assert(face_is_extreme(p, q, r, points)); 
  }
  
  return r; 
}

//return the normal vector of a triangle
point3d normal_vector(triangle3d* t){
  point3d& a = *t->a;
  point3d& b = *t->b;
  point3d& c = *t->c;

  // Calculate the normal vector of the triangle = CB X CA
  point3d BA; // = b - a;
  BA.x = b.x - a.x;
  BA.y = b.y - a.y;
  BA.z = b.z - a.z;
  point3d CA; // = c - a;
  CA.x = c.x - a.x;
  CA.y = c.y - a.y;
  CA.z = c.z - a.z;
  point3d normal; // BA x CA
  normal.x = BA.y * CA.z - BA.z * CA.y;
  normal.y = BA.z * CA.x - BA.x * CA.z;
  normal.z = BA.x * CA.y - BA.y * CA.x;

  return normal;
}

  


// finds and returns a face on the hull3d 
/*
triangle3d find_first_face(vector<point3d>& points) {
  
  edge3d e = find_first_edge_on_hull(points);
  int first_point = e.ia;
  int second_point = e.ib;
    
  //first_point and second_point are both on the hull. Find the third point similarly.
  int third_point = pivot_around_edge(first_point, second_point, points);
  
  printf("%15s", "third point: ");
  print_point(points[third_point], third_point);
  printf("\n"); 

  triangle3d t;
  t.a = &points[first_point];
  t.b = &points[second_point];
  t.c = &points[third_point];

  //if indices are stored 
  t.ia = first_point;
  t.ib = second_point;
  t.ic = third_point;
  
  orient_triangle(&t, points); 
  return t; 
}
*/





/* compute the convex hull of the points */
/*
void gift_wrapping_hull(vector<point3d>& points, vector<triangle3d>& hull) {
  
  hull.clear(); //to be safe
   
  // find a face guaranteed to be on the hull
  triangle3d firstFace = find_first_face(points);
  hull.push_back(firstFace);
  for (int i = 0; i<hull.size(); i++){
    printf("p1:(%d, %d, %d) p2:(%d, %d, %d) p3:(%d, %d, %d)\n", firstFace.a->x, firstFace.a->y, firstFace.a->z, firstFace.b->x, firstFace.b->y, firstFace.b->z, firstFace.c->x, firstFace.c->y, firstFace.c->z);
  }
  
  // find an edge e of a face f that's on the CH, such that the face on the other side of e has not been found 
} 
*/







/*
//////////////////////////////////////////////////////////////////////

                               INCREMENTAL 
//////////////////////////////////////////////////////////////////////
*/

bool point_inside_convex_polygon(point3d p, vector<triangle3d>& hull){
  // if p is inside the polygon it will be behind all the faces (since they should have been oriented facing outwards)
  for (int i = 0; i<hull.size(); i++){
    if (infront_strictly(*(hull[i].a), *(hull[i].b), *(hull[i].c), p)){
      return false;
    }
  }
  return true;
}

// return true if p1 and p2 are identical
bool compare_3dpoints (point3d p1, point3d p2){
  return ((p1.x == p2.x) && (p1.y == p2.y) && (p1.z == p2.z));
}

// returns true if face shares all its edges with those in hull
/*bool shares_all_edges(triangle3d face, vector<triangle3d> hull, int faceIndex){
  printf("entered shares_all_edges\n");
  printf("face: (%d, %d, %d) (%d, %d, %d) (%d, %d, %d)\n", face.a->x, face.a->y, face.a->z, face.b->x, face.b->y, face.b->z, face.c->x, face.c->y, face.c->z);
  for (int i = 0; i<hull.size(); i++){
    printf("hull at start: hull[%d]:(%d, %d, %d) (%d, %d, %d) (%d, %d, %d)\n", i, hull[i].a->x, hull[i].a->y, hull[i].a->z, hull[i].b->x, hull[i].b->y, hull[i].b->z, hull[i].c->x, hull[i].c->y, hull[i].c->z);
  }
  bool edge1_found = false;
  bool edge2_found = false;
  bool edge3_found = false;
  for (int i = 0; i < hull.size(); i++) { // for each face in the hull
    if (i == faceIndex) continue; //  we don't want to consider the face itself when checking if it shares its edges with the hull - it will always be true if we don't do this
    triangle3d cur = hull[i];
    if ((compare_3dpoints(*(cur.a), *(face.a)) && compare_3dpoints(*(cur.b), *(face.b))) || 
        (compare_3dpoints(*(cur.a), *(face.b)) && compare_3dpoints(*(cur.b), *(face.a)))) {
      printf("true: edge1 1 found\n");
      edge1_found = true;
    }
    if ((compare_3dpoints(*(cur.a), *(face.b)) && compare_3dpoints(*(cur.b), *(face.c))) || 
        (compare_3dpoints(*(cur.a), *(face.c)) && compare_3dpoints(*(cur.b), *(face.b)))) {
      printf("true: edge2 1 found\n");
      edge2_found = true;
    }
    if ((compare_3dpoints(*(cur.a), *(face.a)) && compare_3dpoints(*(cur.b), *(face.c))) || 
        (compare_3dpoints(*(cur.a), *(face.c)) && compare_3dpoints(*(cur.b), *(face.a)))) {
      printf("true: edge1 3 found\n");
      edge3_found = true;
    }
    if (edge1_found && edge2_found && edge3_found) { // all edges have been found
      printf("true 1\n");
      return true;
    }
    // Check the next edge
    if ((compare_3dpoints(*(cur.b), *(face.a)) && compare_3dpoints(*(cur.c), *(face.b))) || 
        (compare_3dpoints(*(cur.b), *(face.b)) && compare_3dpoints(*(cur.c), *(face.a)))) {
      printf("true: edge1 2 found\n");
      edge1_found = true;
    }
    if ((compare_3dpoints(*(cur.b), *(face.b)) && compare_3dpoints(*(cur.c), *(face.c))) || 
        (compare_3dpoints(*(cur.b), *(face.c)) && compare_3dpoints(*(cur.c), *(face.b)))) {
      printf("true: edge2 2 found\n");
      edge2_found = true;
    }
    if ((compare_3dpoints(*(cur.b), *(face.a)) && compare_3dpoints(*(cur.c), *(face.c))) || 
        (compare_3dpoints(*(cur.b), *(face.c)) && compare_3dpoints(*(cur.c), *(face.a)))) {
      printf("true: edge3 2 found\n");
      edge3_found = true;
    }
    if (edge1_found && edge2_found && edge3_found) { // all edges have been found
      printf("true 2\n");
      return true;
    }
    // Check the final edge
    if ((compare_3dpoints(*(cur.c), *(face.a)) && compare_3dpoints(*(cur.a), *(face.b))) || 
        (compare_3dpoints(*(cur.c), *(face.b)) && compare_3dpoints(*(cur.a), *(face.a)))) {
      printf("true: edge1 3 found\n");
      edge1_found = true;
    }
    if ((compare_3dpoints(*(cur.c), *(face.b)) && compare_3dpoints(*(cur.a), *(face.c))) || 
        (compare_3dpoints(*(cur.c), *(face.c)) && compare_3dpoints(*(cur.a), *(face.b)))) {
      printf("true: edge2 3 found\n");
      edge2_found = true;
    }
    if ((compare_3dpoints(*(cur.c), *(face.a)) && compare_3dpoints(*(cur.a), *(face.c))) ||
        (compare_3dpoints(*(cur.c), *(face.c)) && compare_3dpoints(*(cur.a), *(face.a)))) {
      printf("true: edge3 3 found\n");
      edge3_found = true;
        }
    if (edge1_found && edge2_found && edge3_found) { // all edges have been found
      printf("true 3\n");
      return true;
    }
  }
  printf("false\n");
  return false;
}*/

// return true if edge is shared with any of the edges of face
bool shares_edge(pair<point3d, point3d> edge, triangle3d face){
  if ((compare_3dpoints(edge.first, *(face.a)) && compare_3dpoints(edge.second, *(face.b))) || 
      (compare_3dpoints(edge.first, *(face.b)) && compare_3dpoints(edge.second, *(face.a))) || 
      (compare_3dpoints(edge.first, *(face.b)) && compare_3dpoints(edge.second, *(face.c))) || 
      (compare_3dpoints(edge.first, *(face.c)) && compare_3dpoints(edge.second, *(face.b))) || 
      (compare_3dpoints(edge.first, *(face.a)) && compare_3dpoints(edge.second, *(face.c))) ||
      (compare_3dpoints(edge.first, *(face.c)) && compare_3dpoints(edge.second, *(face.a)))) {
    return true;
  }
  return false;
}


// return true if edge1 < edge 2
// ordered by x, then y, then z
bool edgesCmp(edge3d edge1, edge3d edge2){
  if (compare_3dpoints(*(edge1.a), *(edge2.a))){
    if (edge1.b->x < edge2.b->x){
      return true;
    }
    else if (edge1.b->x == edge2.b->x){
      if (edge1.b->y < edge2.b->y){
        return true;
      }
      else if (edge1.b->y == edge2.b->y){
        if (edge1.b->z < edge2.b->z){
          return true;
        }
        else{
          return false;
        }
      }
    return false;
  }
  }
  if (edge1.a->x < edge2.a->x){
    return true;
  }
  else if (edge1.a->x == edge2.a->x){
    if (edge1.a->y < edge2.a->y){
      return true;
    }
    else if (edge1.a->y == edge2.a->y){
      if (edge1.a->z < edge2.a->z){
        return true;
      }
      else{
        return false;
      }
  }
  return false;
  }
  return false;
}

// return true if e1 is the same edge as e2
bool sameEdge(edge3d e1, edge3d e2){
  point3d* p11 = e1.a;
  point3d* p21 = e1.b;
  point3d* p12 = e2.a;
  point3d* p22 = e2.b;

  return (p11->x == p12->x && p11->y == p12->y && p11->z == p12->z && p21->x == p22->x && p21->y == p22->y && p21->z == p22->z);
}


void print_hull(vector<triangle3d> hull){
  for (int i = 0; i<hull.size(); i++){
    //printf("new hull: hull[%d]:(%d, %d, %d) hull[%d]:(%d, %d, %d) hull[%d]:(%d, %d, %d)\n", i, hull[i].a->x, hull[i].a->y, hull[i].a->z, i, hull[i].b->x, hull[i].b->y, hull[i].b->z, i, hull[i].c->x, hull[i].c->y, hull[i].c->z);
    printf("(%d, %d, %d) ", hull[i].ia, hull[i].ib, hull[i].ic);
  }
  printf("\n");
}

bool is_extreme_limit(int i, int j, int k, int limit, vector<point3d>& points){
  int nfront=0, nback=0; 
  for (int l =0; l < limit; l++) {
    if (l==i || l==j || l==k)  continue;
    if (coplanar(points[i], points[j], points[k], points[l])) continue; 
    if (infront_strictly(points[i], points[j], points[k], points[l]))
      nfront++; 
    else  nback++; 
    
    if (nfront == 0 || nback  == 0) return false; 
  }
  //if we got here, all on same side 
  return true;
}

void orient_first_face(triangle3d* t, point3d p){
  // make sure that face (0,1,2) faces away from point 3
  point3d* a = t->a;
  point3d* b = t->b;
  point3d* c = t->c;
  if (infront_strictly(*a, *b, *c, p)){
    //printf("orienting triangle\n");
    std::swap(t->b, t->c);
    std::swap(t->ib, t->ic);
  }
}

// orient so that the normal is outside the hull 
// ie all the points in the hull are behind the face or coplanar with it
void orient_triangle(triangle3d* t, vector<triangle3d>& hull) {
  point3d* a = t->a;
  point3d* b = t->b;
  point3d* c = t->c;
  vector<point3d> hullpts;
  for (int i = 0; i<hull.size(); i++){
    hullpts.push_back(*(hull[i].a));
    hullpts.push_back(*(hull[i].b));
    hullpts.push_back(*(hull[i].c));
  }
  //print_triangle(*t, "triangle before orientation", "\n");
  // want to orient t so that all the points in the hull are behind it
  for (int i = 0; i<hullpts.size(); i++){
    if (coplanar(*a, *b, *c, hullpts[i])){
      continue;
    }
    else if (infront_strictly(*a, *b, *c, hullpts[i])){
      //printf("orienting triangle\n");
      std::swap(t->b, t->c);
      std::swap(t->ib, t->ic);
      break;
    }
    else{
      // behind so everything else must also be behind
      break;
    }
  }
  //print_triangle(*t, "triangle after orientation", "\n");
  /*for (int i = 0; i<hullpts.size(); i++){
    // confirm all the points in the hull are either coplanar or behind the face
    //assert(!infront_strictly(*(t->a), *(t->b), *(t->c), hullpts[i]) || coplanar(*(t->a), *(t->b), *(t->c), hullpts[i]));
    long long vol = signed_volume(*(t->a), *(t->b), *(t->c), hullpts[i]);
    if (vol < 0){
      printf("i = %d is in front the face\n", i);
    }
    else if (vol > 0){
      printf("i = %d is behind of the face\n", i);
    }
    else{
      printf("i = %d is coplanar with the face\n", i);
    }
    }*/
  }

/* compute the convex hull of the points */
void incremental_hull(vector<point3d>& points, vector<triangle3d>& hull) {

  for (int i = 0; i<points.size(); i++){
    printf("point[%d]: (%d, %d, %d)\n", i, points[i].x, points[i].y, points[i].z);
  }

  hull.clear(); // to be safe

  //initialise hull to contain faces made up of the first 4 points - no guarantee any of these faces will feature at the end
  triangle3d face1;
  face1.a = &points[0]; face1.ia = 0;
  face1.b = &points[1]; face1.ib = 1;
  face1.c = &points[2]; face1.ic = 2;
  face1.color[0] = (double)rand() / RAND_MAX; // R component
  face1.color[1] = (double)rand() / RAND_MAX; // G component
  face1.color[2] = (double)rand() / RAND_MAX; // B component
  orient_first_face(&face1, points[3]);
  assert(!infront_strictly(*(face1.a), *(face1.b), *(face1.c), points[3]));
  print_triangle(face1, "adding face to hull", "\n");
  hull.push_back(face1);
  
  // once the first face is oriented properly, the remaining faces do not need to be oriented using a function
  // instead we can just flip the shared edge and connect that to the new point

  //take the face containing points (0,2,3)
  //the shared edge is is (0,2) so the correctly oriented face is 

  triangle3d face2;
  face2.a = &points[0]; face2.ia = 0;
  face2.b = &points[2]; face2.ib = 2;
  face2.c = &points[3]; face2.ic = 3;
  face2.color[0] = (double)rand() / RAND_MAX; // R component
  face2.color[1] = (double)rand() / RAND_MAX; // G component
  face2.color[2] = (double)rand() / RAND_MAX; // B component
  orient_triangle(&face2, hull);
  assert(!infront_strictly(*(face2.a), *(face2.b), *(face2.c), points[1]));
  print_triangle(face2, "adding face to hull", "\n");
  hull.push_back(face2);
  //printf("triangle before orientation: (%d, %d, %d)\n", face2.ia, face2.ib, face2.ic);
  //printf("triangle after orientation: (%d, %d, %d)\n", face2.ia, face2.ib, face2.ic);


  triangle3d face3;
  face3.a = &points[1]; face3.ia = 1;
  face3.b = &points[2]; face3.ib = 2;
  face3.c = &points[3]; face3.ic = 3;
  face3.color[0] = (double)rand() / RAND_MAX; // R component
  face3.color[1] = (double)rand() / RAND_MAX; // G component
  face3.color[2] = (double)rand() / RAND_MAX; // B component
  orient_triangle(&face3, hull);
  assert(!infront_strictly(*(face3.a), *(face3.b), *(face3.c), points[0]));
  print_triangle(face3, "adding face to hull", "\n");
  hull.push_back(face3);
  //printf("triangle before orientation: (%d, %d, %d)\n", face3.ia, face3.ib, face3.ic);
  //printf("triangle after orientation: (%d, %d, %d)\n", face3.ia, face3.ib, face3.ic);


  triangle3d face4;
  face4.a = &points[0]; face4.ia = 0;
  face4.b = &points[1]; face4.ib = 1;
  face4.c = &points[3]; face4.ic = 3;
  face4.color[0] = (double)rand() / RAND_MAX; // R component
  face4.color[1] = (double)rand() / RAND_MAX; // G component
  face4.color[2] = (double)rand() / RAND_MAX; // B component
  orient_triangle(&face4, hull);
  assert(!infront_strictly(*(face4.a), *(face4.b), *(face4.c), points[2]));
  print_triangle(face4, "adding face to hull", "\n");
  hull.push_back(face4);
  //printf("triangle before orientation: (%d, %d, %d)\n", face4.ia, face4.ib, face4.ic);
  //printf("triangle after orientation: (%d, %d, %d)\n", face4.ia, face4.ib, face4.ic);


  /*for (int i = 0; i<hull.size(); i++){
    // if the first four starting faces are extreme print true
    if(is_extreme_limit(hull[i].ia, hull[i].ib, hull[i].ic, 4, points)){
      printf("true\n");
    }
  }*/

  for (int i = 0; i<hull.size(); i++){
    print_triangle(hull[i], "current hull", "\n");
  }

  for (int i = 4; i<points.size(); i++){ // for each remaining point p in order 
    // test if p is inside hull
    printf("point %d\n", i);
    if (point_inside_convex_polygon(points[i], hull)){
      printf("point p: (%d, %d, %d) inside hull\n", points[i].x, points[i].y, points[i].z);
      continue; // if p is inside the hull we can ignore it since it definitely does not feature in the hull
    }
    // p is outside of the current hull so we need to add p to the hull and remove all the faces visible from p
    vector<triangle3d> savedFaces;
    for (int j = 0; j<hull.size(); j++){ // for each face f in the hull
      printf("\tprocessing face %d: (%d, %d, %d) \n", j, hull[j].ia, hull[j].ib, hull[j].ic);
      if (infront_strictly(*(hull[j].a), *(hull[j].b), *(hull[j].c), points[i])){ // check if f is visible from p
        // don't include hull[j] in savedFaces
        //printf("visible: erasing (%d, %d, %d), (%d, %d, %d), (%d, %d, %d)\n", hull[j].a->x, hull[j].a->y, hull[j].a->z, hull[j].b->x, hull[j].b->y, hull[j].b->z, hull[j].c->x, hull[j].c->y, hull[j].c->z);
        printf("\t\tvisible so not added to saved faces\n");
        continue; // don't add to savedFaces
      }
      else{
        printf("\t\tnot visible, adding to saved faces\n");
        savedFaces.push_back(hull[j]);
      }
    }

    hull = savedFaces;
    //printf("all faces removed for this current point\n");
    //printf("new hull: \n");
    //print_hull(hull);


    // having removed the visible faces there is now a hole in the hull 
    // figure out boundary of edges around the hole and then add all the new faces - each edge on boundary connected to p

    // add all the edges of the hull (having had visible faces removed) to a vector
    // sort it 
    // all non boundary edges will be in pairs, all boundary edges will be solo, find the solo ones and add them to boundary edges
    vector<edge3d> boundary_edges;
    vector<edge3d> hullEdges;
    for (int i = 0; i<hull.size(); i++){
      edge3d one, two, three;
      one.a = hull[i].a; one.ia = hull[i].ia;
      one.b = hull[i].b; one.ib = hull[i].ib;
      two.a = hull[i].a; two.ia = hull[i].ia;
      two.b = hull[i].c; two.ib = hull[i].ic;
      three.a = hull[i].b; three.ia = hull[i].ib;
      three.b = hull[i].c; three.ib = hull[i].ic;
      hullEdges.push_back(one);
      hullEdges.push_back(two);
      hullEdges.push_back(three);
    }

  /*for (int i = 0; i<hullEdges.size(); i++){
    printf("pre sort and pre orientation: hullEdges[%d]: (%d, %d, %d) (%d, %d, %d)\n", i, hullEdges[i].a->x, hullEdges[i].a->y, hullEdges[i].a->z, hullEdges[i].b->x, hullEdges[i].b->y, hullEdges[i].b->z);
  }*/

    // loop through hullEdges and ensure that they are all oriented in the same way i.e a->b rather than b->a
    // We chose to orient the edges so that the point with the smaller x coordinate comes first 
    for(int j = 0; j<hullEdges.size(); j++){
      point3d* m = hullEdges[j].a;
      point3d* n = hullEdges[j].b;
      if (m->x > n->x){
        hullEdges[j].a = n;
        hullEdges[j].b = m;
      }
    }

  /*for (int i = 0; i<hullEdges.size(); i++){
    printf("pre sort and post orientation: hullEdges[%d]: (%d, %d, %d) (%d, %d, %d)\n", i, hullEdges[i].a->x, hullEdges[i].a->y, hullEdges[i].a->z, hullEdges[i].b->x, hullEdges[i].b->y, hullEdges[i].b->z);
  }*/

    sort(hullEdges.begin(), hullEdges.end(), edgesCmp);

  /*for (int i = 0; i<hullEdges.size(); i++){
    printf("post sort and post orientation: hullEdges[%d]: (%d, %d)\n", i, hullEdges[i].ia, hullEdges[i].ib);
  }*/

  // now we have a vector (hullEdges) containing all the edges in the hull sorted in order and oriented the same way
  // we can now find the boundary edges by identifying the edges that do not come as a pair
  for (int k = 0; k < hullEdges.size(); k++) {
    if (k == 0) { // this needs to be a separate case for 0 and size -1
      if (!sameEdge(hullEdges[k], hullEdges[hullEdges.size() - 1]) && !sameEdge(hullEdges[0], hullEdges[1])) {
        boundary_edges.push_back(hullEdges[k]);
      }
    }
    else if (k == hullEdges.size() - 1){
      if (!sameEdge(hullEdges[k], hullEdges[k-1]) && !sameEdge(hullEdges[k], hullEdges[0])){
        boundary_edges.push_back(hullEdges[k]);
      }
    }
    else if (!sameEdge(hullEdges[k], hullEdges[k+1]) && !sameEdge(hullEdges[k], hullEdges[k-1])) {
      boundary_edges.push_back(hullEdges[k]);
    }
  }


  for (int i = 0; i<boundary_edges.size(); i++){
    printf("boundary_edges[%d]: (%d, %d)\n", i, boundary_edges[i].ia, boundary_edges[i].ib);
  }

  // now we should have boundary edges complete so we can make triangle3d's for each edge to point p (points[i])
  // add new faces to the hull by connecting each boundary edge to the new point
    for (int j = 0; j < boundary_edges.size(); j++) {
      triangle3d new_face;
      new_face.a = boundary_edges[j].a; new_face.ia = boundary_edges[j].ia;
      new_face.b = boundary_edges[j].b; new_face.ib = boundary_edges[j].ib;
      new_face.c = &points[i]; new_face.ic = i;
      new_face.color[0] = (double)rand() / RAND_MAX; // R component
      new_face.color[1] = (double)rand() / RAND_MAX; // G component
      new_face.color[2] = (double)rand() / RAND_MAX; // B component
      //print_triangle(new_face, "adding new face to hull PRE ORIENTATION", "\n");
      orient_triangle(&new_face, hull);
      //print_triangle(new_face, "adding new face to hull POST ORIENTATION", "\n");
      hull.push_back(new_face);
      /*if (is_extreme_limit(new_face.ia, new_face.ib, new_face.ic, i, points)){
        printf("new face added is extreme\n");
      }
      else{
        printf("not extreme\n");
      }*/
      //assert(is_extreme_limit(new_face.ia, new_face.ib, new_face.ic, i, points));
    }
  }

  /*for (int i = 0; i<hull.size(); i++){
    if(!face_is_extreme(hull[i].ia, hull[i].ib, hull[i].ic, points)){
      print_triangle(hull[i], "face", " is not extreme\n");
      //exit(1);
    }
  }*/
}

