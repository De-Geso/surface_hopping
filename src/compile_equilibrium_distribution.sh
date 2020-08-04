# Brad Friesen
# Compile friction.f90
# Updated: 2020-05-11
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

./src/compile_modules.sh

echo 'Building equilibrium_distribution...'
fbuild -Ibin/ -c src/equilibrium_distribution.f90 -o bin/equilibrium_distribution.o
fbuild -llapack -o bin/equilibrium_distribution.out \
	bin/functions.o \
	bin/parameters.o \
	bin/range_finder.o \
	bin/string_utilities.o \
	bin/equilibrium_distribution.o
echo 'equilibrium_distribution built'
