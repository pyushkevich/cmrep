// ********************************************************
//
// Source file for DrawFillInside
//
// ********************************************************



void 
draw_fill_inside_image(unsigned char *image, int *dims)
  // Fill the 'inside' part of an image (closed objects that do not touch the
  // image edges)
{
  int x, y, z;
  int size_plane = dims[0]*dims[1];
  int size_line = dims[0];

  int foreground = 255;
  int label = 1;
  int background = 0;

  // fill image edges -> won't work if object touches the edge !!
  for (z = 0; z < dims[2]; z++) {
    for (y = 0; y < dims[1]; y++) {
      if (y == 0 || z == 0) { // draw whole plane
	for (x = 0; x < dims[0] ; x++) 
	  image[x +  size_line * y + size_plane * z] = label;
      } else { // draw edges of x
	image[0 +  size_line * y + size_plane * z] = label;
	image[size_line - 1 +  size_line * y + size_plane * z] = label;
      }
    }
  }
  
  // forward propagation
  for (z = 1; z < dims[2]-1; z++) {
    for (y = 1; y < dims[1]-1; y++) {
      for (x = 1; x < dims[0]-1; x++) {
	int index = x +  size_line * y + size_plane * z;
	if (image[index] == background &&    // check past neighborhood
	    (image[index - 1] == label || 
	     //image[index - size_line + 1 ] == label || 
	     image[index + size_line] == label || 
	     //image[index - size_line - 1 ] == label == label ||
	     //image[index - size_plane - size_line - 1 ] == label ||
	     //image[index - size_plane - size_line] == label ||
	     //image[index - size_plane - size_line + 1 ] == label ||
	     //image[index - size_plane - 1] == label ||
	     image[index - size_plane] == label 
	     //image[index - size_plane + 1] == label || 
	     //image[index - size_plane + size_line - 1 ] == label ||
	     //image[index - size_plane + size_line] == label ||
	     //image[index - size_plane + size_line + 1 ] == label 
	     )) {
	  image[index] = label;
	}
      }
    }
  }

  // backward propagation
  for (z = dims[2]-1; z > 0; z--) {
    for (y = dims[1]-1; y > 0; y--) {
      for (x = dims[0]-1; x > 0; x--) {
	int index = x +  size_line * y + size_plane * z;
	if (image[index] == background &&    // check past neighborhood
	    (image[index + 1] == label || 
	     //image[index + size_line + 1 ] == label || 
	     image[index + size_line] == label || 
	     //image[index + size_line - 1 ] == label ||
	     //image[index + size_plane - size_line - 1 ] == label ||
	     //image[index + size_plane - size_line] == label ||
	     //image[index + size_plane - size_line + 1 ] == label ||
	     //image[index + size_plane - 1] == label ||
	     image[index + size_plane] == label 
	     //image[index + size_plane + 1] == label || 
	     //image[index + size_plane + size_line - 1 ] == label ||
	     //image[index + size_plane + size_line] == label ||
	     //image[index + size_plane + size_line + 1 ] == label 
	     )) {
	  image[index] = label;
	}
      }
    }
  }

  // reassign labels
  for (int i = 0; i < dims[2]*dims[1]*dims[0]; i++) {
    if (image[i] == label) image[i] = background;
    else if (image[i] == background) image[i] = foreground;
  }
}
