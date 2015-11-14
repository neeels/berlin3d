#!/usr/bin/env python
import sys, os, os.path, time
from filelock import *

# For benefit of speed, this parses for a very specific formatting of XML found
# in these particular datasets.

DEST_DIR = 'gml-sections'

def section_of(point):
  return tuple([(int)(c / 1000.) for c in point])

def section_str(point):
  return '%04d_%04d_%04d' % section_of(point)

def parse_floats(floats_str):
  tokens = floats_str.split()
  return [float(t) for t in tokens]

def midpoint(floats_lists):
  n = float(len(floats_lists))
  return tuple([sum(v)/n for v in zip(*floats_lists)])



def store(src_dir, mid, content):
  content = content.replace('imageURI>appearance',
                            'imageURI>%sappearance' % src_dir)
  fn = '%s/%s' % (DEST_DIR, section_str(mid))
  with FileLock(fn):
    f = open(fn, 'a')
    f.write(content)
    f.close()


def handle_src_dir(src_dir):
  if not src_dir.endswith('/'):
    src_dir += '/'

  bezirk = os.path.basename(src_dir)
  gml_fname = src_dir + 'citygml.gml'

  if not os.path.exists(gml_fname):
    print "Invalid source:", gml_fname
    exit(1)

  if not os.path.exists(DEST_DIR):
    os.makedirs(DEST_DIR)

  print gml_fname, '-->', DEST_DIR

  member_str = None
  member_lower_corner = None
  member_lower_upper = None

  for line in open(gml_fname):
    if line.endswith('\r\n'):
      line = line[:-2] + '\n'

    if line.startswith(' <cityObjectMember>'):
      member_str = [ '<cityObjectMember bezirk="%s">' % bezirk ]
      continue

    elif line.startswith(' </cityObjectMember>'):
      member_str.append(line)

      store(src_dir,
            midpoint((parse_floats(member_lower_corner),
                      parse_floats(member_upper_corner))),
            ''.join(member_str))

      member_str = None
      member_lower_corner = None
      member_lower_upper = None
      exit(0)
      continue


    # spaces match only the bldg:Building/gml:boundedBy/gml:Envelope/gml:{lower,upper}Corner
    if line.startswith('     <gml:lowerCorner>'):
      member_lower_corner = line[22:-19]
    elif line.startswith('     <gml:upperCorner>'):
      member_upper_corner = line[22:-19]

    if member_str is not None:
      member_str.append(line)
      

for src_dir in sys.argv[1:]:
  handle_src_dir(src_dir)

#vim: expandtab nocin ai
