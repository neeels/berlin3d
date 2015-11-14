#!/bin/sh
do_all=0
if [ "$1" = "-c" ]; then
  shift
  do_all=1
fi

target() {
  echo "jpgs/$(basename "$1").jpg"
}

one_ecw2jpg() {
  rm out.rgb
  mkdir -p jpgs
  set -e
  test -f "$1"
  make ecw2raw
	LD_LIBRARY_PATH=./lib ./ecw2raw src-ecw/dop20_404_5810.ecw
  test -f out.rgb
  convert -size 8192x8192 -depth 8 out.rgb -quality 60 "$(target "$1")"
}

for f in $@; do
  if [ "$do_all" = "1" ]; then
    if [ -f "$(target "$f")" ]; then
      continue
    fi
  fi
  one_ecw2jpg "$f"
done

