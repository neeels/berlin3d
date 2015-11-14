#!/usr/bin/env python
import sys, os, os.path

# For benefit of speed, this parses for a very specific formatting of XML found
# in these particular datasets.

DEST_DIR = 'gml-sections'

def section_of(point):
  return tuple([(int)(c / 1000.) for c in point])

def section_str(point):
  return '%04d_%04d_%04d' % section_of(point)

def parse_floats(floats_str):
  tokens = floats_str.split();
  return [float(t) for t in tokens]

def midpoint(floats_lists):
  n = float(len(floats_lists))
  return tuple([sum(v)/n for v in zip(*floats_lists)])

def store(mid, tags):
  fn = '%s/%s' % (DEST_DIR, section_str(mid))
  f = open(fn, 'a')
  f.write(tags)
  f.close()

name = sys.argv[1]
src_dir = 'src-names/' + name + '/'
gml_fname = src_dir + 'citygml.gml';

if os.path.exists(DEST_DIR):
  print 'Adding to %s' % DEST_DIR
else:
  print 'Creating %s' % DEST_DIR
  os.makedirs(DEST_DIR)

member_str = None
member_lower_corner = None
member_lower_upper = None

for line in open(gml_fname):
  if line.startswith(' <cityObjectMember>'):
    member_str = [ line ]

  # spaces match only the bldg:Building/gml:boundedBy/gml:Envelope/gml:{lower,upper}Corner
  elif line.startswith('     <gml:lowerCorner>'):
    member_lower_corner = line[22:-20]
  elif line.startswith('     <gml:upperCorner>'):
    member_upper_corner = line[22:-20]

  elif line.startswith(' </cityObjectMember>'):
    member_str.append(line)
    store(midpoint((parse_floats(member_lower_corner),
                    parse_floats(member_upper_corner))),
          ''.join(member_str))
    member_str = None
    member_lower_corner = None
    member_lower_upper = None
    

#vim: expandtab nocin ai
