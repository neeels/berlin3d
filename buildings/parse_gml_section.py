#!/usr/bin/env python

LIMIT_BUILDINGS = None
DEFINED_MID = (387690.658211, 5801906.009784, 94.15)

import sys, os.path
from xml.sax import saxutils, handler, make_parser

DEST_DIR = 'cpp-sections'

src_section = sys.argv[1]
name = os.path.basename(src_section)

def cross(a, b):
      return [a[1] * b[2] - a[2] * b[1],
              a[2] * b[0] - a[0] * b[2],
              a[0] * b[1] - a[1] * b[0]]

def sub(a, b):
  return [a[0]-b[0],
          a[1]-b[1],
          a[2]-b[2]]

def parse_floats(floats_str, n=3):
  tokens = floats_str.split();
  if n and len(tokens) != n:
    print "This string does not have %d numbers:\n%r" % (n, floats_str)
    raise RuntimeError()
  return [float(t) for t in tokens]

def mid_and_size(floats_lists):
  minv = list(floats_lists[0])
  maxv = list(floats_lists[0])
  for floats in floats_lists:
    for i in range(len(floats)):
      minv[i] = min(minv[i], floats[i])
      maxv[i] = max(maxv[i], floats[i])

  mid = [0.] * len(floats_lists[0])
  for i in range(len(mid)):
    mid[i] = (maxv[i] + minv[i])/2.;

  size = [(maxv[i] - minv[i]) for i in range(len(minv))]

  return mid, size

def mid_and_offsets(floats_lists):
  mid, size = mid_and_size(floats_lists)
  return mid, offsets(mid, floats_lists)

def offsets(mid, floats_lists):
  ofs = []
  for floats in floats_lists:
    ofs.append([(floats[i] - mid[i]) for i in range(len(floats))])

  return ofs


class Building:
  id_str = None
  lower_corner = None
  upper_corner = None
  pos = None
  size = None
  walls = None

  def parse_corners(s):
    s.pos, s.size = mid_and_size((s.lower_corner, s.upper_corner))

  def valid(s):
    s.parse_corners()
    valid = ((s.lower_corner) and (s.upper_corner)
             and(len(s.lower_corner) == 3)
             and (len(s.upper_corner) == 3)
             and (s.size[0] * s.size[1] * s.size[2] > 0.))

    #if not valid:
    #     print "Invalid: %r %r %r" % (s.id_str, s.lower_corner, s.upper_corner)
    return valid

class Done(Exception):
	pass

class CityGml(handler.ContentHandler):

    def __init__(self, out = sys.stdout):
        handler.ContentHandler.__init__(self)
        self._out = out

    def startDocument(self):
        self.buildings = []
        self.building = None
        self.depth = 0
        self.state = 0
        self.building_depth = -1
        self.tex_coords = {}
        self.tex_id = None
        self.ring_id = None
        
    def startElement(self, name, attrs):
      self.depth += 1
      self.contents = ''
      if name == 'bldg:Building':
        self.building_depth = self.depth
        self.building = Building()
        self.building.id_str = attrs.get('gml:id')
        self.building.walls = []
        self.state = 1

      elif self.state == 1 and self.depth == (self.building_depth+1) and name == 'gml:boundedBy':
        self.state = 2

      elif self.state == 2 and name == 'gml:lowerCorner':
        self.state = 3

      elif self.state == 2 and name == 'gml:upperCorner':
        self.state = 4

      elif name == 'app:textureCoordinates':
        self.tex_id = attrs.get('ring')

      elif name == 'gml:LinearRing':
        self.ring_id = attrs.get('gml:id')


    def endElement(self, name):
      if self.building and name == 'gml:posList':
        self.building.walls.append((self.ring_id, parse_floats(self.contents, 0)))

      elif self.state == 1 and name == 'bldg:Building':
        if self.building and self.building.valid():
          self.buildings.append(self.building)
          if LIMIT_BUILDINGS:
            if len(self.buildings) >= LIMIT_BUILDINGS:
              raise Done()
        self.building = None
        self.state = 0

      elif self.state == 2 and name == 'gml:boundedBy':
        self.state = 1

      elif self.state == 3 and name == 'gml:lowerCorner':
        self.building.lower_corner = parse_floats(self.contents, 3)
        self.state = 2

      elif self.state == 4 and name == 'gml:upperCorner':
        self.building.upper_corner = parse_floats(self.contents, 3)
        self.state = 2

      elif name == 'app:textureCoordinates':
        tex_entry = (self.image_uri, parse_floats(self.contents, 0))
        self.tex_coords[self.tex_id] = tex_entry
        self.tex_id = None

      elif name == 'gml:linearRing':
        self.ring_id = None

      elif name == 'app:imageURI':
        self.image_uri = self.contents

      self.depth -= 1

    def characters(self, content):
      self.contents += content

gml = CityGml()

parser = make_parser()
parser.setContentHandler(gml)
try:
	parser.parse(src_section)
except Done:
  pass


h = open(name + '.h', 'w')
h.write('''#pragma once

#include <city.h>
extern const DwellingData {0};

extern const BuildingData {0}_buildings[];
extern const double {0}_points[][3];
extern const double {0}_normals[][3];
extern const double {0}_tex_coords[][2];
extern const int {0}_wall_indices[];
extern const char *{0}_images[];
'''.format(name)
  )
h.close()

cpp = open(name + '.cpp', 'w')

sys.stdout = cpp

print '/* Data generated from citygml: %d buildings */' % len(gml.buildings)
print '''
#include "%s.h"
''' % name

print '''
const DwellingData %s = {
  .name = "%s",
''' % (name, name)

points = []
tex_coords = []
wall_indices = []
normals = []
wall_count = 0;
images = []

for b in gml.buildings:
  b.walls_start_idx = len(wall_indices)
  b.normals_start_idx = len(normals)
  for ring_id, r in b.walls:
    wall_points = []
    for i in range(0, len(r), 3):
      p = r[i:i+3]
      wall_points.append(p)

    tex_points = []
    wall_tex = gml.tex_coords.get('#' + ring_id)
    if ((wall_tex is not None)
        and (len(wall_tex[1])/2 == len(wall_points))):
      img, pts = wall_tex
      wall_indices.append(len(images))
      images.append(img)
      wall_indices.append(len(tex_coords))
      for i in range(0, len(pts), 2):
        tex_coords.append(pts[i:i+2])
    else:
      wall_indices.append(-1)
      wall_indices.append(-1)

    for p in wall_points:
      wall_indices.append(len(points))
      points.append(p)

    wall_indices.append(-2)

    normals.append(cross(sub(wall_points[1], wall_points[0]),
                         sub(wall_points[2], wall_points[1])))


    wall_count += 1

if DEFINED_MID is None:
  mid, points = mid_and_offsets(points)
else:
  mid = DEFINED_MID
  points = offsets(mid, points)

points_min = list(points[0])
points_max = list(points[0])
for p in points:
  points_min[0] = min(points_min[0], p[0])
  points_min[1] = min(points_min[1], p[1])
  points_min[2] = min(points_min[2], p[2])
  points_max[0] = max(points_max[0], p[0])
  points_max[1] = max(points_max[1], p[1])
  points_max[2] = max(points_max[2], p[2])


ofs = offsets(mid, [b.pos for b in gml.buildings])
for i in range(len(gml.buildings)):
  gml.buildings[i].pos = ofs[i]

print '  .mid = {%f, %f, %f},' % tuple(mid)
print '  .points_min = {%f, %f, %f},' % tuple(points_min)
print '  .points_max = {%f, %f, %f},' % tuple(points_max)
print '  .n_buildings = %d,' % len(gml.buildings)
print '  .buildings = %s_buildings,' % name
print '  .n_points = %d,' % len(points)
print '  .points = %s_points,' % name
print '  .n_walls = %d,' % wall_count
print '  .n_wall_indices = %d,' % len(wall_indices)
print '  .wall_indices = %s_wall_indices,' % name
print '  .n_normals = %d,' % len(normals)
print '  .normals = %s_normals,' % name
print '  .n_tex_coords = %d,' % len(tex_coords)
print '  .tex_coords = %s_tex_coords,' % name
print '  .n_images = %d,' % len(images)
print '  .images = %s_images,' % name
print '};'

print 'const BuildingData %s_buildings[] = {' % name
print ',\n'.join([
  ('  { '
   + ',\n    '.join((
          #'.id = "%s"' % b.id_str,
          '.pos = {%f, %f, %f}' % tuple(b.pos),
          '.box_size = {%f, %f, %f}' % tuple(b.size),
          '.walls_start_idx = %d' % b.walls_start_idx,
          '.normals_start_idx = %d' % b.normals_start_idx,
          '.n_walls = %d' % len(b.walls),
          )
          )
   + '  }')
  for b in gml.buildings])

print '''
};
'''

point_strs=['{%f, %f, %f}'%tuple(p) for p in points]

print 'const double %s_points[][3] = {' % name
print ',\n'.join(point_strs)
print '};'

print 'const int %s_wall_indices[] = {' % name
print ',\n'.join([str(i) for i in wall_indices])
print '};'

normals_strs=['{%f, %f, %f}'%tuple(p) for p in normals]
print 'const double %s_normals[][3] = {' % name
print ',\n'.join(normals_strs)
print '};'

tex_coords_strs=['{%f, %f}'%tuple(p) for p in tex_coords]
print 'const double %s_tex_coords[][2] = {' % name
print ',\n'.join(tex_coords_strs)
print '};'

images_strs=['"%s%s"'%(src_dir, s) for s in images]
print 'const char *%s_images[] = {' % name
print ',\n'.join(images_strs)
print '};'

sys.stdout = sys.__stdout__
#vim: expandtab nocin ai
