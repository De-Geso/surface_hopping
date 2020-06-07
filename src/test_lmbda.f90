! test_lmbda.f90
! Brad Friesen
! Look at some simple functions, such as the potentials, and probability
! distributions
! Updated: 2020-05-11

program test_lmbda
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
real lmbda
integer i

call init()

open (unit=10, file='running.dat')
do i = 1, n
	lmbda = lmbda1 + dlmbda*(i-1.)
	write (10,*) lmbda, force_variance(lmbda)
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

end program
