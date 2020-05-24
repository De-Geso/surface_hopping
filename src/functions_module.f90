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
! ENERGIES
!=======================================================================

function V0 (x, lmbda)
! System when an electron has been transferred from metal to molecule.
real V0, x, lmbda
V0 = 0.5*k0*x**2+(G0+lmbda)
end function

function V1 (x)
! System when the molecule is uncharged.
real V1, x
V1 = 0.5*k1*(x-g)**2
end function

function pmf (x, lmbda)
! Potential of mean force
real pmf, x, lmbda
pmf = -1./beta * log(exp(-beta*V0(x, lmbda))+exp(-beta*V1(x)))
end function

function free_energy (lmbda)
! Free energy
real free_energy, lmbda
free_energy = -1./beta * log(sqrt(2.*PI/beta)*(1./sqrt(k1)+exp(-beta*(G0+lmbda))/sqrt(k0)))
end function

!=======================================================================
! PROBABILITIES
!=======================================================================

function prob0 (x, lmbda)
! Probability on V0 potential
real prob0, x, lmbda
prob0 = exp(-beta*V0(x,lmbda))/		&
	(sqrt(2.*pi/(beta*k0))*exp(-beta*(G0+lmbda)) + sqrt(2.*pi/(beta*k1)))
end function

function prob1 (x, lmbda)
! Probability on V1 potential
real prob1, x, lmbda
prob1 = exp(-beta*V1(x))/		&
	(sqrt(2.*pi/(beta*k0))*exp(-beta*(G0+lmbda)) + sqrt(2.*pi/(beta*k1)))
end function

function prob_eq (x, lmbda)
! Probability distribution
real prob_eq, x, lmbda
prob_eq = exp(beta*(free_energy(lmbda)-pmf(x,lmbda)))
end function

function cum_prob_eq (x, lmbda)
! Cumulative probability distribution
real cum_prob_eq, x, lmbda
cum_prob_eq = (sqrt(k1)*(erf(sqrt(beta*k0/2.)*x)+1.) &
	+ sqrt(k0)*exp(beta*(G0+lmbda))*(erf(sqrt(beta*k1/2.)*(x-g))+1.)) &
	/ (2.*(sqrt(k0)*exp(beta*(G0+lmbda))+sqrt(k1)))
end function

function dl_cum_prob_eq (x, lmbda)
! Derivative of cumulative probability distribution with respect to
! lambda.
real dl_cum_prob_eq, x, lmbda
dl_cum_prob_eq = -beta*sqrt(k0*k1)*exp(beta*lmbda) &
	* (erf(sqrt(beta*k0/2.)*x)-erf(sqrt(beta*k1/2.)*(x-g))) &
	/ (2.*(sqrt(k0)*exp(beta*lmbda)+sqrt(k1))**2)
end function


function fermi(e)
real fermi, e
fermi = 1./(1.+exp(beta*e))
end function
!=======================================================================
! INTEGRANDS
!=======================================================================

function int_x_prob_eq (x)
real int_x_prob_eq, x
int_x_prob_eq = prob_eq(x, lmbda0)
end function

function int_x_friction (x)
real int_x_friction, x
int_x_friction = 1./D * (dl_cum_prob_eq(x, lmbda0)**2 &
	/prob_eq(x, lmbda0))
end function

end module
