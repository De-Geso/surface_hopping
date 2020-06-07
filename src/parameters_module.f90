! parameters.f90
! Brad Friesen
! Compile with: gfortran -O3 -fdefault-real-8 -c parameters.f90
! Updated 2020-05-11

module parameters
implicit none

! Physical constants
real, parameter :: Pi=4.D0*DATAN(1.D0)

! System parameters that we don't expect to change
! Spatial offset
real, parameter :: g = 1.0
! Diffusion coefficient
real, parameter :: D = 1.0
! Tolerance for range finder
real, parameter :: tol = 0.001

! Control file parameters
! Force constants
real k0, k1
! Test values
real x0, lmbda0
! Protocol values
real tau
! Range where stuff usually happens
real x_min, x_max, lmbda_min, lmbda_max
! Energy offset
real G0
! Thermodynamic beta
real beta

end module
