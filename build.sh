#!/bin/bash
# Build the SemiSim 3D core. -fno-finite-math-only keeps the NaN guards
# (x!=x tests in the solver) alive under -ffast-math.
set -e
cd "$(dirname "$0")"
if [ -f part1.c ]; then cat part1.c part2a.c part2b.c part2c.c part3.c part3g.c part4.c > semiconductor_3d.c; fi
gcc -O3 -march=native -ffast-math -fno-finite-math-only -fopenmp -shared -fPIC \
    -o semiconductor_core.so semiconductor_3d.c -lm
echo "built semiconductor_core.so"
