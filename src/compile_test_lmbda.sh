# Brad Friesen
# Compile test_lmbda.f90
# Updated: 2020-06-03
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

./src/compile_modules.sh

echo 'Building test_lmbda...'
fbuild -Ibin/ -c src/test_lmbda.f90 -o bin/test_lmbda.o
fbuild -o bin/test_lmbda.out \
	bin/integration.o \
	bin/interpolation.o \
	bin/functions.o \
	bin/parameters.o \
	bin/string_utilities.o \
	bin/range_finder.o \
	bin/root_finding.o \
	bin/test_lmbda.o
echo 'test_lmbda built'
