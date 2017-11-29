#!/usr/bin/env python

LIMIT_BUILDINGS = None

import sys, os.path
from xml.sax import saxutils, handler, make_parser

name = sys.argv[1]
gml_fname = name + '-data/citygml.gml';


class Done(Exception):
	pass

class CityGml(handler.ContentHandler):

    def __init__(self, out = sys.stdout):
        handler.ContentHandler.__init__(self)
        self._out = out

    def startDocument(self):
			self.tags = {}
        
    def startElement(self, name, attrs):
			self.tags[name] = self.tags.get(name, 0) + 1

gml = CityGml()

parser = make_parser()
parser.setContentHandler(gml)
try:
	parser.parse(gml_fname)
except Done:
  pass

for name, nr in sorted(gml.tags.items()):
  print '%6d %s' % (nr, name)
#vim: expandtab nocin ai
