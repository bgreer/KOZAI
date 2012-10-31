
orbit : orbit.F90 bsint.F90 physics.F90
	gfortran -c physics.F90
	gfortran -c bsint.F90
	gfortran -fopenmp orbit.F90 bsint.F90 physics.F90 -o orbit
