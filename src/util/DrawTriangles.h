// ********************************************************
//
// Header file for DrawTriangles.cc
//
// ********************************************************

#ifndef _DRAW_TRIANGLES_H_
#define _DRAW_TRIANGLES_H_

void drawBinaryTrianglesFilled(unsigned char *image, int *dim, 
			       double ** vertex_table,
			       int num_triangles);
//Scan converts a closed object whose surface is described by triangles, the vertex_table
// always groups 3 points to triangles, each point is an array of 3 doubles

void drawBinaryPolygonsFilled(unsigned char *image, int *dim,
			      double** vertex_polygon_table,
			      int num_polygons, int *num_points_per_polygon);
//Scan converts a closed object described by a set of polygons
// the vertex_polygon_table stores all points of the num_polygons polygons
// in  num_points_per_polygon[i] the number of points of the i-th polygon has to be stored

void drawBinaryPolygonsFilled(unsigned char *image, int *dim,
			      double** vertex_polygon_table,
			      int num_polygons, int *num_points_per_polygon,
			      double ** &vertex_triangle_table, int &numTriangles);
//Scan converts a closed object described by a set of polygons
// the vertex_polygon_table stores all points of the num_polygons polygons
// in  num_points_per_polygon[i] the number of points of the i-th polygon has
// to be stored. Returns in vertex_triangle_table the triangle table



void drawBinaryTrianglesSheetFilled(unsigned char *image, int *dim, 
				    double ** vertex_table,
				    int num_triangles);
//Scan converts a manifold object described by a set of triangles
// the coordinates of vertex_table are in units of voxels with 
// image(0,0,0) as an origin
// vertex_table is a table of num_triangles*3 Pointers to 3 doubles (a point)
     
void drawBinaryPolygonsSheetFilled(unsigned char *image, int *dim,
				   double** vertex_polygon_table,
				   int num_polygons, int *num_points_per_polygon);
//Scan converts a manifold object described by a set of polygons
// the vertex_polygon_table stores all points of the num_polygons polygons
// in  num_points_per_polygon[i] the number of points of the i-th polygon has to be stored

void 
drawBinaryPolygonsSheetFilled(unsigned char *image, int *dim,
			      double** vertex_polygon_table,
			      int num_polygons, int *num_points_per_polygon,
			      double ** &vertex_triangle_table, int &numTriangles);
//Scan converts a manifold object described by a set of polygons
// the vertex_polygon_table stores all points of the num_polygons polygons
// in  num_points_per_polygon[i] the number of points of the i-th polygon has
// to be stored. Returns in vertex_triangle_table the triangle table
#endif
