# Brad Friesen
# Compile friction.f90
# Updated: 2020-05-11
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

./src/compile_modules.sh

echo 'Building fpe...'
fbuild -Ibin/ -c src/fpe.f90 -o bin/fpe.o
fbuild -llapack -o bin/fpe.out \
	bin/integration.o \
	bin/interpolation.o \
	bin/linear_algebra.o \
	bin/functions.o \
	bin/parameters.o \
	bin/range_finder.o \
	bin/string_utilities.o \
	bin/fpe.o
echo 'fpe built'
