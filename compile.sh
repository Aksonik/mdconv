#!/bin/bash
#gcc -c amber.c
#gcc -c charmm.c
#gcc -c datastream.c
#gcc -c format.c
#gcc -c trajdata.c
#gcc -static -O3 -o mdconv mdconv.c *.o -lm -lquadmath

g++ -c amber.c
g++ -c charmm.c
g++ -c datastream.c
g++ -c format.c
g++ -c trajdata.c
g++ -c mdconv.c
g++ -c pdb.c
g++ -o mdconv *.o
#g++ -static -O3 -o mdconv mdconv.c *.o -lm -lquadmath
