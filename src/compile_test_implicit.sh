# Compile test_implicit.f90
# Updated: 2020-05-11
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

./src/compile_modules.sh

echo 'Building test_implicit...'
fbuild -Ibin/ -c src/test_implicit.f90 -o bin/test_implicit.o
fbuild -llapack -o bin/test_implicit.out		\
	bin/linear_algebra.o		\
	bin/test_implicit.o
echo 'test_implicit built'
