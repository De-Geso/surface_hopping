# Brad Friesen
# Compile all the modules
# Updated: 2020-05-11
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

echo 'Building interpolation...'
fbuild -Jbin/ -o bin/interpolation.o -c src/modules/interpolation.f90
echo 'interpolation built'

echo 'Building integration...'
fbuild -Jbin/ -o bin/integration.o -c src/modules/integration.f90
echo 'integration built'

echo 'Building linear_algebra...'
fbuild -Jbin/ -o bin/linear_algebra.o -c src/modules/linear_algebra.f90
echo 'linear_algebra built'

echo 'Building parameters...'
fbuild -Jbin/ -o bin/parameters.o -c src/parameters_module.f90
echo 'parameters built'

echo 'Building string_utilities'
fbuild -Jbin/ -o bin/string_utilities.o -c src/string_utilities_module.f90
echo 'string_utilities built'

echo 'Building functions...'
fbuild -Jbin/ -o bin/functions.o -c src/functions_module.f90
echo 'functions built'
