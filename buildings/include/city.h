#pragma once

#include "pt.h"

typedef unsigned int uint;

struct BuildingData {
//  const char *const id;
  Pt pos;
  Pt box_size;
  uint walls_start_idx;
  uint normals_start_idx;
  uint n_walls;
};

struct DwellingData {
	const char *name;
  Pt mid;
  Pt points_min;
  Pt points_max;
  uint n_buildings;
	const BuildingData *buildings;
  uint n_points;
	const double (*points)[3];
  uint n_walls;

	/* wall_indices is of format:

		 <images-idx|-1> <tex_coords-idx|-1> <point-idx> <point-idx> [...] -2  [ <images-idx... ]

		 e.g. for a wall with no texture:
		   -1 -1 0 1 2 3 4 -2
		 or with texture:
		   2 3 0 1 2 3 4 -2
  */
  uint n_wall_indices;
	const int *wall_indices;
  uint n_normals;
	const double (*normals)[3];
  uint n_tex_coords;
  const double (*tex_coords)[2];
  uint n_images;
  const char **const images;
};

