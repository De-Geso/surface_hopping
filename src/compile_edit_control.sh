# Brad Friesen
# Compile path.f90
# Updated: 2020-06-04
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

./src/compile_modules.sh

echo 'Building edit_control...'
fbuild -Ibin/ -c src/edit_control.f90 -o bin/edit_control.o
fbuild -o bin/edit_control.out \
	bin/parameters.o \
	bin/string_utilities.o \
	bin/edit_control.o
echo 'edit_control built'
