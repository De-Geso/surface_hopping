! force_variance_path.f90
! Brad Friesen
! Approximate work, calculated using the force variance
! Updated: 2020-05-11

program force_variance_work
use integration
use interpolation
use functions
use parameters
use range_finder
use string_utilities
implicit none


character(:), allocatable :: datafile, outfile, tempfile
integer, parameter :: n = 1 + 100
real, dimension(:), allocatable :: x, y, y2
real lmbda, alpha, optimal, naive, avg, sqrtavg
integer nlines
integer i


call init()

! Find average values of force variance for approximate ratio
call qromb(int_avg, lmbda_min, lmbda_max, avg)
avg = avg/(lmbda_max-lmbda_min)
call qromb(int_sqrtavg, lmbda_min, lmbda_max, sqrtavg)
sqrtavg = sqrtavg/(lmbda_max-lmbda_min)

! Find work
call qromb(int_naive, 0., tau, naive)
call qromb(int_optimal, 0., tau, optimal)

! Record results
open (unit=10, file=outfile, position='append')
write (10,*) tau, optimal, naive, avg/sqrtavg**2
close(10)


contains


subroutine init()
integer i
real yp1, ypn
character(:), allocatable :: datafile

! Read control file variables
call read_control()
call find_lrange()
call find_xrange()

! Find length of optimal path file
datafile = make_filename('data/path/', 'path')
nlines = filelength(datafile)
allocate(x(nlines))
allocate(y(nlines))
allocate(y2(nlines))

! Populate arrays with optimal path information
open(1, file=datafile)
do i = 1, nlines
	read(1, *) x(i), y(i)
end do
close(1)

! Spline the optimal path
yp1 = (y(2)-y(1))/(x(2)-x(1))
ypn = (y(n)-y(n-1))/(x(n)-x(n-1))
call spline (x, y, nlines, yp1, ypn, y2)

! Integrate alpha
call qromb(int_alpha, lmbda_min, lmbda_max, alpha)
alpha = alpha/tau

! Make filenames
tempfile = make_filename('data/work/', 'temp')
outfile = make_filename('data/work/', 'force_variance_work')
end subroutine



!=======================================================================
! INTEGRANDS
!=======================================================================
function int_avg (lmbda)
real int_avg, lmbda
int_avg = force_variance(lmbda)
end function

function int_sqrtavg (lmbda)
real int_sqrtavg, lmbda
int_sqrtavg = sqrt(force_variance(lmbda))
end function

function int_naive (t)
real int_naive, t, v, lmbda
v = (lmbda_max-lmbda_min)/tau
lmbda = lmbda_min + v*t
int_naive = v*v*force_variance(lmbda)
write (40,*) t/tau, lmbda
end function

function int_optimal (t)
real int_optimal, t, v, lmbda
call splint(x, y, y2, nlines, t/tau, lmbda)
v = alpha/sqrt(force_variance(lmbda))
int_optimal = v*v*force_variance(lmbda)
end function

function int_alpha (lmbda)
real int_alpha, lmbda
int_alpha = sqrt(force_variance(lmbda))
end function

end program
 
