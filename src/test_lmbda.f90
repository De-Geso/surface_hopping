! test_lmbda.f90
! Brad Friesen
! Look at some functions, such as the potentials, and probability
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
character(:), allocatable :: test_file, ef_file, hopping_file, marcus_file
real lmbda, test(3), ef(3), hopping(3), marcus(3), fwd, rev
integer i


call init()

open (unit=100, file=test_file)
open (unit=10, file=ef_file)
open (unit=11, file=hopping_file)
open (unit=12, file=marcus_file)
do i = 1, n
	lmbda = lmbda1 + dlmbda*(i-1.)


	! Surface Hopping Marcus theory case
	! (relaxation time from Marcus theory)
	! Friction
	call qromb(int_test_friction, x_min, x_max, test(1))
	! Force Variance (is analytic)
	test(2) = force_variance(lmbda)
	! Relaxation Time
	test(3) = test_relax(lmbda)
	write(100,*) lmbda, test(:)
		
	! Electron Friction
	! Friction
	call qromb(int_ef_friction, x_min, x_max, ef(1))
	! Force Variance
	! First moment
	call qromb(int_ef_conj_1, x_min, x_max, ef(2))
	! Second moment
	call qromb(int_ef_conj_2, x_min, x_max, ef(3))
	ef(2) = ef(3) - ef(2)*ef(2)
	! Relaxation Time
	ef(3) = ef(1)/ef(2)
	write(10,*) lmbda, ef(:)

	! Surface Hopping hopping case
	! (relaxation time from hopping rates only)
	! Friction
	call qromb(int_hopping_friction, x_min, x_max, hopping(1))
	! Force Variance (is analytic)
	hopping(2) = force_variance(lmbda)
	! Relaxation Time
	hopping(3) = hopping(1)/hopping(2)
	write(11,*) lmbda, hopping(:)

	! Surface Hopping Marcus theory case
	! (relaxation time from Marcus theory)
	! Friction
	call qromb(int_marcus_friction, x_min, x_max, marcus(1))
	! Force Variance (is analytic)
	marcus(2) = force_variance(lmbda)
	! Relaxation Time
	marcus(3) = marcus_relax(lmbda)
	write(12,*) lmbda, marcus(:)
end do



contains



subroutine init()
call read_control()
call find_lrange()
call find_xrange()
test_file = make_filename ('data/lmbda/', 'lmbda_test')
ef_file = make_filename ('data/lmbda/', 'lmbda_ef')
hopping_file = make_filename ('data/lmbda/', 'lmbda_hopping', gam=gam)
marcus_file = make_filename ('data/lmbda/', 'lmbda_marcus', gam=gam)
end subroutine



!=======================================================================
! TESTING
!=======================================================================
real function int_test_friction (x)
real x, fwd, rev
! Friction integrand
int_test_friction = 1./(4.*cosh(beta/2.*(V(1,x,lmbda)-V(0,x,lmbda)))**2) &
	*test_relax(lmbda)*(prob_eq(0,x,lmbda)+prob_eq(1,x,lmbda))
write (1,*) lmbda, 1.0/(fwd+rev)
end function


real function test_relax (lmbda)
real lmbda, fwd, rev
! Relaxation time
! Forward electron transfer rate
call qromb(int_test_fwd, x_min, x_max, fwd)
! Reverse electron transfer rate
call qromb(int_test_rev, x_min, x_max, rev)
test_relax = 1.0/(fwd+rev)
end function


real function int_test_fwd (x)
real x
! Forward rate from equilibrium average rates
int_test_fwd = gam*dt*(1.-fermi(V(1, x, lmbda)-V(0, x, lmbda)))*prob_eq(0, x, lmbda)
end function


real function int_test_rev (x)
real x
! Forward rate from equilibrium average rates
int_test_rev = gam*dt*(fermi(V(1, x, lmbda)-V(0, x, lmbda)))*prob_eq(1, x, lmbda)
end function


!=======================================================================
! SURFACE HOPPING (Marcus)
!=======================================================================
real function int_marcus_friction (x)
real x, fwd, rev
! Integrand for the SH friction, using Marcus theory
! Whole shabang
int_marcus_friction = 1./(4.*cosh(beta/2.*(V(1,x,lmbda)-V(0,x,lmbda)))**2) &
	*marcus_relax(lmbda)*(prob_eq(0,x,lmbda)+prob_eq(1,x,lmbda))
write (1,*) lmbda, 1.0/(fwd+rev)
end function


real function marcus_relax (lmbda)
real lmbda, fwd, rev
! Relaxation time, using Marcus Theory
! Forward electron transfer rate
fwd = 2.*pi*(gam*dt)**2/sqrt(4.*pi/beta*(V(1,0.,lmbda)-V(1,1.,lmbda))) &
	* exp(-(V(1,0.,lmbda)-V(1,1.,lmbda)-lmbda)**2/(4./beta*(V(1,0.,lmbda)-V(1,1.,lmbda))))
! Reverse electron transfer rate
rev = 2.*pi*(gam*dt)**2/sqrt(4.*pi/beta*(V(0,1.,lmbda)-V(0,0.,lmbda))) &
	* exp(-(V(0,1.,lmbda)-V(0,0.,lmbda)+lmbda)**2/(4./beta*(V(0,1.,lmbda)-V(0,0.,lmbda))))
marcus_relax = 1.0/(fwd+rev)
end function



!=======================================================================
! SURFACE HOPPING (hopping)
!=======================================================================
real function int_hopping_friction (x)
real x
! Integrand for the SH friction, using Miranda's derivation
int_hopping_friction = beta/(4.*gam*dt*cosh(beta/2.*(V(1,x,lmbda)-V(0,x,lmbda)))**2) &
	* (prob_eq(0,x,lmbda)+prob_eq(1,x,lmbda))
end function


!=======================================================================
! ELECTRON FRICTION
!=======================================================================
real function int_ef_friction (x)
real x
! Integrand for the EF friction, using Zulkowski
int_ef_friction = dl_cum_prob_eq(x,lmbda)**2/(prob_eq(0,x,lmbda)+prob_eq(1,x,lmbda))
end function


real function int_ef_conj_1 (x)
real x
int_ef_conj_1 = ef_conjugate_force(x, lmbda)*(prob_eq(0,x,lmbda)+prob_eq(1,x,lmbda))
end function


real function int_ef_conj_2 (x)
real x
int_ef_conj_2 = ef_conjugate_force(x, lmbda)**2*(prob_eq(0,x,lmbda)+prob_eq(1,x,lmbda))
end function


real function ef_conjugate_force (x, lmbda)
real x, lmbda
! Force conjugate to the control parameter lmbda
ef_conjugate_force = exp(-beta*V(0,x,lmbda))/(exp(-beta*V(0,x,lmbda))+exp(-beta*V(1,x,lmbda)))
end function


real function dl_cum_prob_eq (x, lmbda)
! Derivative of cumulative probability distribution with respect to
! lambda.
real x, lmbda
dl_cum_prob_eq = -beta*sqrt(k0*k1)*exp(beta*lmbda) &
	* (erf(sqrt(beta*k0/2.)*x)-erf(sqrt(beta*k1/2.)*(x-g))) &
	/ (2.*(sqrt(k0)*exp(beta*lmbda)+sqrt(k1))**2)
end function



!real function marcus_relaxation(x)
!real f, x, fwd, rev
!! Forward rate
!fwd = 2.*pi*gam*gam &
!	/ sqrt(4.*pi*(1./beta)*(V1(0.)-V1(1.))) &
!	* exp(-(V1(0.)-V1(1.)-lmbda)**2/(4./beta*(V1(0.)-V1(1.))))
!! Reverse rate
!rev = 2.*pi*gam*gam &
!	/ sqrt(4.*pi*(1./beta)*(V0(0.,lmbda)-V0(1.,lmbda))) &
!	* exp(-(V0(0.,lmbda)-V0(1.,lmbda)-lmbda)**2/(4./beta*(V0(0.,lmbda)-V0(1.,lmbda))))
!marcus_relax = 1./(fwd+rev)*(prob0(x, lmbda)+prob1(x, lmbda))
!end function



!pure function pmf (x, lmbda)
!intent (in) x, lmbda
!real pmf, x, lmbda
!! Potential of mean force of the function. Free energy averaged over
!! both potentials.
!pmf = -1.0/beta * log(exp(-beta*V0(x,lmbda)) + exp(-beta*V1(x)))
!end function



!function pmf_force_variance (lmbda)
!real pmf_force_variance, lmbda, s(2)
!call qromb(int_fvar1, x_min, x_max, s(1))
!call qromb(int_fvar2, x_min, x_max, s(2))
!pmf_force_variance = s(2) - s(1)*s(1)
!end function




end program
