#!/bin/sh
#
# Convert all given ecw files to jpg.
#
#  -c  skip if target jpg already exists (default: overwrite)
skip_existing=0
if [ "$1" = "-c" ]; then
  shift
  skip_existing=1
fi

target() {
  echo "jpgs/$(basename "$1").jpg"
}

one_ecw2jpg() {
  mkdir -p jpgs
  set -e
  test -f "$1"
	src="$1"
  dest="$(target "$1")"
	LD_LIBRARY_PATH=./lib ./ecw2jpg "$src" "$dest"
  echo "wrote $dest"
}

for f in $@; do
  if [ "$skip_existing" = "1" ]; then
    if [ -f "$(target "$f")" ]; then
      continue
    fi
  fi
  one_ecw2jpg "$f"
done

