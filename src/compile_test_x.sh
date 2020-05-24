# Brad Friesen
# Compile test_x.f90
# Updated: 2020-05-11
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

./src/compile_modules.sh

echo 'Building test_x...'
fbuild -Ibin/ -c src/test_x.f90 -o bin/test_x.o
fbuild -o bin/test_x.out \
	bin/integration.o bin/interpolation.o \
	bin/functions.o bin/parameters.o bin/string_utilities.o bin/test_x.o
echo 'test_x built'
