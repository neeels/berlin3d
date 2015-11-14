#!/bin/sh
mkdir -p cpp-sections
(

echo "all: \\"

cd gml-sections
for section in *; do
  echo " ${section}.cpp \\"
done
echo

echo "%.cpp: ../gml-sections/%"
echo "\t./parse_gml_section.py \$@"

) > cpp-sections/Makefile

make -C cpp-sections

