! range_finder.f90
! author: Brad Friesen
! updated: 2019-11-09
! purpose: Find an appropriate lambda range and xrange for given k's and G0
! Compile: gfortran -O3 -fdefault-real-8 range_finder.f90
! Run: run_range_finder.sh

module range_finder
use root_finding
use parameters
use string_utilities
use functions
implicit none

contains

subroutine find_lrange ()
lmbda_max = -log(sqrt(k0/k1)*(tol/(1.0-tol))) / beta
lmbda_min = -log(sqrt(k0/k1)*(1.0/tol-1.0)) / beta

call edit_ctrl ('lmbda_min', lmbda_min)
call edit_ctrl ('lmbda_max', lmbda_max)
end subroutine find_lrange



subroutine find_xrange ()
real z0, z1, a, b

z0 = sqrt(2.0*pi/(beta*k0))
z1 = sqrt(2.0*pi/(beta*k1))

a = sqrt(-2.0/(beta*k0)*log(tol*z0))
b = sqrt(-2.0/(beta*k1)*log(tol*z1))

x_min = min(-a, -b+g)
x_max = max(a, b+g)

call edit_ctrl ('x_min', x_min)
call edit_ctrl ('x_max', x_max)
end subroutine find_xrange

end module range_finder
