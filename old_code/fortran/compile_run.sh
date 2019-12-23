#!/bin/bash

rm fpres/*.dat
rm err/*.dat
gfortran -c fancy_timestep.f90
gfortran -o runfp runfp.f90 fancy_timestep.o -llapack -lfftw3
./runfp
