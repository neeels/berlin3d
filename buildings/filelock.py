import os, time
import errno

class FileLock:
  def __init__(s, path):
    s.path = path
    s.dirname = os.path.dirname(path)
    s.basename = os.path.basename(path)
    s.lock = path + '.locked'
    s.acquired = 0

  def __enter__(s):
    while s.acquired < 1:
      try:
        os.mkdir(s.lock)
        s.acquired += 1
        return
      except OSError as e:
        if e.errno != errno.EEXIST:
          raise
        print 'waiting for lock:', s.lock
        time.sleep(.1)
        pass

  def __exit__(s, *ign_args, **ign_kwargs):
    if s.acquired == 1:
      os.rmdir(s.lock)
      s.acquired = 0
    else:
      s.acquired -= 1

  def __del__(s):
    if s.acquired > 0:
      os.rmdir(s.lock)
      s.acquired = 0
