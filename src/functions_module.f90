! functions.f90
! Brad Friesen
! Functions, such as potentials and probability distributions
! Compile with: gfortran -O3 -fdefault-real-8 -c functions.f90
! Updated: 2020-05-10

module functions
use parameters
use string_utilities
implicit none

contains

!=======================================================================
! ENERGETICS
!=======================================================================
! Energies, forces, etc.

pure function V (alpha, x, lmbda)
intent (in) alpha, x, lmbda
! Potential of the system. alpha=0 for V0, alpha=1 for V1
real V, x, lmbda
integer alpha
V = (0.5*k0*x**2+(G0+lmbda))*(1-alpha) &
	+ (0.5*k1*(x-g)**2)*alpha
end function



!pure function V0 (x, lmbda)
!intent (in) x, lmbda
!! System when an electron has been transferred from metal to molecule.
!real V0, x, lmbda
!V0 = 0.5*k0*x**2+(G0+lmbda)
!end function



!pure function V1 (x)
!intent (in) x
!! System when the molecule is uncharged.
!real V1, x
!V1 = 0.5*k1*(x-g)**2
!end function



!pure function E (x, lmbda)
!intent(in) x, lmbda
!! Energy gap between the states
!real E, x, lmbda
!E = V1(x) - V0(x, lmbda)
!end function



pure function force_variance (lmbda)
intent (in) lmbda
real force_variance, lmbda
! Equilibrium average conjugate force variance of the system at a fixed
! lambda. Conjugate force is trivially -1 when on V0, and 0 when on V1
real a
a = 1.0/(1.0 + sqrt(k0/k1)*exp(beta*(lmbda+G0)))
force_variance = a - a*a
end function



!=======================================================================
! PROBABILITIES
!=======================================================================

pure function free_energy (lmbda)
intent(in) lmbda
! Helmholtz free energy of the system at a particular lambda
real free_energy, lmbda
free_energy = -1./beta*log(partition_function(lmbda))
end function


pure function partition_function (lmbda)
intent(in) lmbda
real partition_function, lmbda
! Partition function
partition_function = sqrt(2.*pi/(beta*k0))*exp(-beta*(G0+lmbda)) + sqrt(2.*pi/(beta*k1))
end function



pure function prob_eq (alpha, x, lmbda)
intent (in) alpha, x, lmbda
! Probability of a particular state of the system.
real prob_eq, x, lmbda
integer alpha
prob_eq = exp(-beta*V(alpha,x,lmbda))/partition_function(lmbda)
end function



!function prob0 (x, lmbda)
!! Probability on V0 potential
!real prob0, x, lmbda
!prob0 = exp(-beta*V0(x,lmbda))/		&
!	(sqrt(2.*pi/(beta*k0))*exp(-beta*(G0+lmbda)) + sqrt(2.*pi/(beta*k1)))
!end function




!function prob1 (x, lmbda)
!! Probability on V1 potential
!real prob1, x, lmbda
!prob1 = exp(-beta*V1(x))/		&
!	(sqrt(2.*pi/(beta*k0))*exp(-beta*(G0+lmbda)) + sqrt(2.*pi/(beta*k1)))
!end function


function cum_prob_eq (x, lmbda)
! Cumulative probability distribution
real cum_prob_eq, x, lmbda
cum_prob_eq = (sqrt(k1)*(erf(sqrt(beta*k0/2.)*x)+1.) &
	+ sqrt(k0)*exp(beta*(G0+lmbda))*(erf(sqrt(beta*k1/2.)*(x-g))+1.)) &
	/ (2.*(sqrt(k0)*exp(beta*(G0+lmbda))+sqrt(k1)))
end function



function fermi(e)
real fermi, e
fermi = 1./(1.+exp(beta*e))
end function

!=======================================================================
! "SINGLE VARIABLE"
!=======================================================================
! Multivariable stuff set as single variable stuff. i.e. we only change
! one variable at a time. For integrands, or searching in 1 dimension.

end module
