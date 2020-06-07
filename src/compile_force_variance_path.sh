# Brad Friesen
# Compile path.f90
# Updated: 2020-05-11
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

./src/compile_modules.sh

echo 'Building force_variance_path...'
fbuild -Ibin/ -c src/force_variance_path.f90 -o bin/force_variance_path.o
fbuild -o bin/force_variance_path.out \
	bin/integration.o \
	bin/interpolation.o \
	bin/functions.o \
	bin/parameters.o \
	bin/range_finder.o \
	bin/string_utilities.o \
	bin/force_variance_path.o
echo 'force_variance_path built'
