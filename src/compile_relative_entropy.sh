# Brad Friesen
# Compile friction.f90
# Updated: 2020-05-11
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

./src/compile_modules.sh

echo 'Building relative_entropy...'
fbuild -Ibin/ -c src/relative_entropy.f90 -o bin/relative_entropy.o
fbuild -llapack -o bin/relative_entropy.out \
	bin/parameters.o \
	bin/string_utilities.o \
	bin/relative_entropy.o
echo 'relative_entropy built'
