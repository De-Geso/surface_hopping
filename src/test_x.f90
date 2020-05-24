! test_x.f90
! Brad Friesen
! Look at some simple functions, such as the potentials, and probability
! distribution.
! Updated: 2020-05-11

program test_x
use functions
use parameters
use string_utilities
use integration
implicit none

real x, xstep, xmax, xmin, sl, s, sr
integer i, n

parameter (n=1000, xmax=g+1., xmin=-1., xstep=(xmax-xmin)/(n-1.))


call init()

do i = 1, n
	x = xmin + xstep*(i-1.)
	write (10,*) x, prob_eq(x, lmbda0), prob0(x,lmbda0), prob1(x,lmbda0)
end do

call qromo(int_x_prob_eq, -1.e30, -10., sl, midinf)
call qromo(int_x_prob_eq, 10., 1e30, sr, midinf)
call qromo(int_x_prob_eq, -10., 10., s, midpnt)
write (*, *) sl+s+sr


contains

subroutine init()
call read_control()
end subroutine

end program
