! force_variance_path.f90
! Brad Friesen
! Approximate work, calculated using the force variance
! Updated: 2020-05-11

program force_variance_work
use functions
use parameters
use range_finder
use string_utilities
implicit none


character(:), allocatable :: outfile
integer, parameter :: n = 1000
real, parameter :: xx_min = -2.0
real, parameter :: xx_max = 3.0
real, parameter :: dx = (xx_max-xx_min)/(n-1.0)

real, dimension(n, 4) :: potential, distribution
real, dimension(n) :: x
real, dimension(2) :: z
integer i


! Read control file variables
call read_control()
call find_lrange()
call find_xrange()

outfile = make_filename('data/distribution/', 'dist_eq')

do i = 1, n
	x(i) = xx_min + dx*(i-1.)
	potential(i, 1) = V0(x(i),lmbda_min)
	potential(i, 2) = V1(x(i))
	potential(i, 3) = V0(x(i),lmbda_max)
	potential(i, 4) = V1(x(i))
end do

do i = 1, n
	z(1) = z(1) + exp(-beta*potential(i,1)) + exp(-beta*potential(i,2))
	z(2) = z(2) + exp(-beta*potential(i,3)) + exp(-beta*potential(i,4))
end do

do i = 1, n
	distribution(i, 1) = exp(-beta*V0(x(i), lmbda_min))/z(1)
	distribution(i, 2) = exp(-beta*V1(x(i)))/z(1)
	distribution(i, 3) = exp(-beta*V0(x(i), lmbda_max))/z(2)
	distribution(i, 4) = exp(-beta*V1(x(i)))/z(2)
end do

open (unit=10, file=outfile)
do i = 1,n
	write (10,*) x(i), distribution(i,1), distribution(i,2), distribution(i,3), distribution(i,4)
end do
close (10)

end program
