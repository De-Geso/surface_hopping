# author: Brad Friesen
# updated: 2020-06-04
# purpose: automate several runs of data
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real- $@
}

# program name that we want to run
program='test_x'

./src/compile_$program.sh
./src/compile_edit_control.sh

k0=(4.0 16.0 64.0)
k1=(4.0 16.0 64.0)
# tau=(0.01 0.1 0.5 1.0 2.0 10.0 50.0 200.0)
# gam=(0.1 1.0 10.0 10000.0)

for i in ${!k0[@]}; do
	echo "Editing control file..."
	./bin/edit_control.out 'k0' ${k0[$i]}

	for j in ${!k1[@]}; do
		echo "Editing control file..."
		./bin/edit_control.out 'k1' ${k1[$j]}
		
		#for k in ${!tau[@]}; do
			#echo "Editing control file..."
			#./bin/edit_control.out 'tau' ${tau[$k]}

			#for l in ${!gam[@]}; do
				#echo "Editing control file..."
				#./bin/edit_control.out 'gam' ${gam[$l]}
                
                echo "Running" $program.f90 '...'
                time ./bin/$program.out
			#done
		#done
	done
done
