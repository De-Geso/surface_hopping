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
real lmbda, s(3), friction, fvar, pmf_fric, pmf_fvar
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
	pmf_fric = pmf_friction(lmbda)
	pmf_fvar = pmf_force_variance(lmbda)
	write (10,*) lmbda, friction, fvar, friction/fvar, pmf_fric, pmf_fvar, pmf_fric/pmf_fvar
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


function pmf_friction (lmbda)
real pmf_friction, lmbda, s(3)
call qromb(int_pmf_friction, x_min, x_max, pmf_friction)
!call qromo(int_pmf_friction, -1.e30, x_min, s(1), midinf)
!call qromo(int_pmf_friction, x_min, x_max, s(2), midpnt)
!call qromo(int_pmf_friction, x_max, 1.e30, s(3), midinf)
!pmf_friction = sum(s)
end function


function int_pmf_friction (x)
real int_pmf_friction, x
int_pmf_friction = dl_cum_prob_eq(x, lmbda)**2/(prob0(x, lmbda)+prob1(x, lmbda))
end function


pure function pmf (x, lmbda)
intent (in) x, lmbda
real pmf, x, lmbda
! Potential of mean force of the function. Free energy averaged over
! both potentials.
pmf = -1.0/beta * log(exp(-beta*V0(x,lmbda)) + exp(-beta*V1(x)))
end function


pure function pmf_conjugate_force (x, lmbda)
intent(in) x, lmbda
real pmf_conjugate_force, x, lmbda
! Force conjugate to the control parameter lmbda
pmf_conjugate_force = exp(-beta*V0(x,lmbda))/ &
	(exp(-beta*V1(x))+exp(-beta*V0(x,lmbda)))
end function


function pmf_force_variance (lmbda)
real pmf_force_variance, lmbda, s(2)
call qromb(int_fvar1, x_min, x_max, s(1))
call qromb(int_fvar2, x_min, x_max, s(2))
!call qromo(int_pmf_force_variance, -1.e30, x_min, s(1), midinf)
!call qromo(int_pmf_force_variance, x_min, x_max, s(2), midpnt)
!call qromo(int_pmf_force_variance, x_max, 1.e30, s(3), midinf)
pmf_force_variance = s(2) - s(1)*s(1)
end function


function int_fvar1 (x)
real int_fvar1, x
int_fvar1 = pmf_conjugate_force(x, lmbda) &
	*(prob0(x, lmbda)+prob1(x, lmbda))
end function


function int_fvar2 (x)
real int_fvar2, x
int_fvar2 = pmf_conjugate_force(x, lmbda)**2 &
	*(prob0(x, lmbda)+prob1(x, lmbda))
end function

end program
