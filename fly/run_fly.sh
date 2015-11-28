#!/bin/sh
./fly -g $(xrandr | grep "\*" | sed 's/^[ \t]*\([0-9x]*\)[ \t]*.*/\1/')
