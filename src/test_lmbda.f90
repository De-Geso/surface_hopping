! test_lmbda.f90
! Brad Friesen
! Look at some simple functions, such as the potentials, and probability
! distributions
! Updated: 2020-05-11

program test_lmbda
use integration
use functions
use parameters
use range_finder
use string_utilities
implicit none

integer, parameter :: n = 1 + 1000
real, parameter :: lmbda1 = -10.
real, parameter :: lmbda2 = 10.
real, parameter :: dlmbda = (lmbda2-lmbda1)/(n-1.)
character(:), allocatable :: outfile
real lmbda, s(3), friction, fvar
integer i

call init()

open (unit=10, file='running.dat')
do i = 1, n
	lmbda = lmbda1 + dlmbda*(i-1.)
	call qromo(f, -1.e30, x_min, s(1), midinf)
	call qromo(f, x_min, x_max, s(2), midpnt)
	call qromo(f, x_max, 1.e30, s(3), midinf)
	friction = sum(s)
	fvar = force_variance(lmbda)
	write (10,*) lmbda, friction, fvar, friction/fvar
end do
close (unit=10)
call rename('running.dat', outfile)



contains



subroutine init()
call read_control()
call find_lrange()
call find_xrange()
outfile = make_filename ('data/lmbda/', 'lmbda')
end subroutine


function f(x)
real f, x
f = 1./cosh(beta*E(x, lmbda)/2.)**2 * (exp(-beta*V0(x,lmbda))+exp(-beta*V1(x)))/Z(lmbda)
end function

end program
