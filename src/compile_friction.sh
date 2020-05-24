# Brad Friesen
# Compile friction.f90
# Updated: 2020-05-11
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

./src/compile_modules.sh

echo 'Building friction...'
fbuild -Ibin/ -c src/friction.f90 -o bin/friction.o
fbuild -o bin/friction.out \
	bin/integration.o bin/interpolation.o \
	bin/functions.o bin/parameters.o bin/string_utilities.o bin/friction.o
echo 'friction built'
