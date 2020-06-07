! force_variance_path.f90
! Brad Friesen
! Optimal path, calculated using the force variance
! Updated: 2020-05-11

program force_variance_path
use integration
use functions
use parameters
use range_finder
use string_utilities
implicit none

integer, parameter :: n = 100000
character(:), allocatable :: outfile, tempfile
real lmbda, alpha, t, dt
integer i


call init()

open (unit=10, file=tempfile)
write (10,*) t/tau, lmbda
do i = 1,n-1
	t = t+dt
	lmbda = lmbda + dt*alpha/sqrt(force_variance(lmbda))
	write (10,*) t/tau, lmbda
end do
! Cheat on the very last point
write (10,*) 1., lmbda_max
close (10)
call rename(tempfile, outfile)


contains


subroutine init()
call read_control()
call find_lrange()
call find_xrange()

lmbda = lmbda_min

call qromb (int_alpha, lmbda_min, lmbda_max, alpha)
alpha = alpha/tau

t = 0.
dt = tau/n

tempfile = make_filename ('data/path/', 'temp')
outfile = make_filename ('data/path/', 'path')
end subroutine


function int_alpha (lmbda)
real int_alpha, lmbda
int_alpha = sqrt(force_variance(lmbda))
end function

end program
