# Brad Friesen
# Compile all the modules
# Updated: 2020-05-11
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

# Numerical recipes modules
fbuild -Jbin/ -o bin/interpolation.o -c src/modules/interpolation.f90
fbuild -Jbin/ -o bin/integration.o -c src/modules/integration.f90
fbuild -Jbin/ -o bin/linear_algebra.o -c src/modules/linear_algebra.f90
fbuild -Jbin/ -o bin/root_finding.o -c src/modules/root_finding.f90

# Other modules
fbuild -Jbin/ -o bin/functions.o -c src/functions_module.f90
fbuild -Jbin/ -o bin/parameters.o -c src/parameters_module.f90
fbuild -Jbin/ -o bin/range_finder.o -c src/range_finder_module.f90
fbuild -Jbin/ -o bin/string_utilities.o -c src/string_utilities_module.f90
