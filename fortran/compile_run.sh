#!/bin/bash

rm fpres/*.dat
gfortran -c solvefp.f90 fancy_timestep.f90
gfortran -o runfp runfp.f90 solvefp.o fancy_timestep.o -llapack -lfftw3
./runfp
