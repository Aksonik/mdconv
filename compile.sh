#!/bin/bash
gcc -c amber.c
gcc -c charmm.c
gcc -c datastream.c
gcc -c format.c
gcc -c trajdata.c
gcc -static -O3 -o mdconv mdconv.c *.o -lm -lquadmath
