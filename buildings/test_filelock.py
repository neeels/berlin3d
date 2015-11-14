#!/usr/bin/env python

import sys

from filelock import *


with FileLock('xx'):
  print 'acquired'
  sys.stdin.readline()

print 'released'

# vim: expandtab nocin ai
