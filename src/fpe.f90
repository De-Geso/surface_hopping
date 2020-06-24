! fpe.f90
! Brad Friesen
! Solve the Fokker Planck equation
! Updated: 2020-05-13

program fpe
use interpolation
use functions
use parameters
use range_finder
use string_utilities
implicit none

! Theta = 1.0 for implicit, theta = 0.5 for Crank-Nicolson
real, parameter :: theta = 0.5
! Grid limits we are solving on
real, parameter :: xx_min = -2.0
real, parameter :: xx_max = 3.0
integer, parameter :: n = 1000
real, parameter :: dt = 0.0001
real, parameter :: dx = (xx_max-xx_min)/(n-1.0)
real, parameter :: eps = 1.e-8
integer, parameter :: check_step = int(1.0/dt)
character(:), allocatable :: outfile, tempfile
! Spline variables
integer nlines
real, dimension(:), allocatable :: x, y, y2
! FPE grids
real, dimension(n) :: positions
real, dimension(n, 2) :: prob, p_now, p_last, p_eq, p_last_ref, &
	pot_at_pos, drift, diff
real z, lmbda, work
integer i, j



call init()
if (gam*dt .gt. 1.) stop 'Hopping constant greater than 1!'
lmbda = lmbda_min
call init_optimal()

! Populate arrays
do i = 1, n
	pot_at_pos(i, 1) = V0(positions(i),lmbda)
	pot_at_pos(i, 2) = V1(positions(i))
end do

! Calculate the partition function
do i = 1, n
do j = 1, 2
	z = z + exp(-beta*pot_at_pos(i, j))
end do
end do

! Calculate the equilibrium distribution
prob = exp(-beta*pot_at_pos)/z
! prob(:,1) = 1./(2.*n)
! prob(:,2) = 1./(2.*n)

p_now = prob

call steady_state()
call evolve()
call dump()
write (*,*) work



contains



subroutine init()
! Initialize variables
integer i

call read_control()
call find_lrange()
call find_xrange()

prob = 0.0
p_now = 0.0
p_last = 0.0
p_last_ref = 0.0
pot_at_pos = 0.0
diff = D
z = 0.0
work = 0.0

outfile = make_filename('data/fpe/', 'fpe', tau=tau)
tempfile = make_filename('data/fpe/', 'temp', tau=tau)

do i = 1, n
	positions(i) = xx_min + dx*(i-1.)
	drift(i,1) = drift0(positions(i))
	drift(i,2) = drift1(positions(i))
end do
end subroutine



subroutine init_optimal()
! Get optimal path data and spline it
integer i
real yp1, ypn
character(:), allocatable :: datafile

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

end subroutine



subroutine steady_state()
real tot_var_dist
integer continue_condition, step_counter
integer i, j

tot_var_dist = 0.0
continue_condition = 1
step_counter = 0

do while (continue_condition .eq. 1)
	! Save previous distribution
	p_last = p_now
	! Advance probability one timestep on a potential
	call update_prob()
	! Slosh probability between potentials
	call update_transfer(lmbda)
	
	if (step_counter .ge. check_step) then
		! check if the distribution changed from last check
		tot_var_dist = 0.5*sum(abs(p_last_ref-p_now))
		write (*,*) tot_var_dist, sum(p_now) - 1.

		! Check normalization
		call check()

		! End loop if we're at equilibrium
		if (tot_var_dist .lt. eps) then
			continue_condition = 0
			! Save current distribution as equilibrium
			p_eq = p_now
		else
			tot_var_dist = 0.0
			step_counter = 0
			p_last_ref = p_now
		end if
	end if
	step_counter = step_counter + 1
end do
end subroutine


subroutine evolve()
real t, lmbda_last
integer nt, i, j

nt = int(tau/dt)

do i = 1, nt
	! Save last lambda value for calculating the work
	lmbda_last = lmbda
	! Naive straight line path
	lmbda = lmbda_min + (lmbda_max-lmbda_min)/tau*(i*dt)
	! Save previous distribution
	p_last = p_now
	! Advance probability one timestep on a potential
	call update_prob()
	! Slosh probability between potentials
	call update_transfer(lmbda)
	! Calculate work
	do j = 1, n
		work = work + &
			(V0(positions(j),lmbda)-V0(positions(j),lmbda_last))*p_now(j,1)
	end do
	! work = work-(free_energy(lmbda)-free_energy(lmbda_last))
	write (100, *) i, work

	call check()
end do
end subroutine



subroutine update_prob()
! Update probability by solving the Fokker-Planck with explicit,
! implicit, or crank-nicholson schemes
real dl(n-1,2), d(n,2), du(n-1,2), b(n,2)
real ll(n-1,2), dd(n,2), uu(n-1,2), bb(n,2)
integer info(2)
integer i

! Construct the tridiagonal matrix, all unknowns on the left
do i = 1, n-1
	du(i,:) = -theta*dt*(diff(i+1,:)/(dx*dx) - drift(i+1,:)/(2.0*dx))
	dl(i,:) = -theta*dt*(diff(i,:)/(dx*dx) + drift(i,:)/(2.0*dx))
end do
d = 1.0 - theta*dt*(-2.0*diff/(dx*dx))

! Construct the solution b, all knowns on the right
do i = 2, n-1
	b(i,:) = p_last(i,:) + dt*(1.0-theta)*( &
		-(drift(i+1,:)*p_last(i+1,:) - drift(i-1,:)*p_last(i-1,:))/(2.0*dx) &
		+(diff(i+1,:)*p_last(i+1,:)-2.0*diff(i,:)*p_last(i,:)+diff(i-1,:)*p_last(i-1,:))/(dx*dx))
end do

! von Neumann boundary conditions
! Left
d(1,:) = 1.0 - theta*dt*(drift(1,:)/(2.0*dx) + (-diff(1,:)-dx*drift(1,:))/(dx*dx))
du(1,:) = -theta*dt*(-drift(2,:)/(2.0*dx) + diff(2,:)/(dx*dx))
b(1,:) = p_last(1,:) + (1.0-theta)*dt*( &
	-(drift(2,:)*p_last(2,:) - drift(1,:)*p_last(1,:))/(2.0*dx) &
	+(diff(2,:)*p_last(2,:) - (diff(1,:)+dx*drift(1,:))*p_last(1,:))/(dx*dx))	
! Right
d(n,:) = 1.0 - theta*dt*(-drift(n,:)/(2.0*dx) + (-diff(n,:)+dx*drift(n,:))/(dx*dx))
dl(n-1,:) = -theta*dt*(drift(n-1,:)/(2.0*dx) + diff(n-1,:)/(dx*dx))
b(n,:) = p_last(n,:) + (1.0-theta)*dt*( &
	-(drift(n,:)*p_last(n,:) - drift(n-1,:)*p_last(n-1,:))/(2.0*dx) &
	+((-diff(n,:)+dx*drift(n,:))*p_last(n,:) + diff(n-1,:)*p_last(n-1,:))/(dx*dx))

! copies to feed to improve if necessary
bb = b
uu = du
dd = d
ll = dl

do i = 1,2
	call dgtsv(n, 1, dl(:,i), d(:,i), du(:,i), b(:,i), n, info(i))
	if (info(i) .NE. 0) stop info(i)
	! p_now(:,i) = improve(n, ll(:,i), dd(:,i), uu(:,i), bb(:,i), b(:,i))
end do
p_now = b
end subroutine



function improve(n, dl, d, du, b, x)
real improve(n), dl(n-1), d(n), du(n-1), b(n), x(n)
integer n
real r(n)
integer i, info
double precision sdp

sdp = -b(1)
sdp = sdp + d(1)*x(1) + du(1)*x(2)
r(1) = sdp
do i = 2, n-1
	sdp = -b(i)
	sdp = sdp + dl(i-1)*x(i-1) + d(i)*x(i) + du(i)*x(i+1)
	r(i) = sdp
end do
sdp = -b(n)
sdp = sdp + dl(n-1)*x(n-1) + d(i)*x(i)
r(n) = sdp

call dgtsv(n, 1, dl, d, du, r, n, info)
if (info .NE. 0) stop info

do i = 1, n
	improve(i) = x(i)-r(i)
end do
end function



subroutine update_transfer(lmbda)
real de, lmbda, trans01, trans10
integer i

do i = 1, n
	de = V1(positions(i))-V0(positions(i), lmbda)
	trans01 = gam*dt*fermi(de)*p_now(i,1)
	trans10 = gam*dt*(1.-fermi(de))*p_now(i,2)
	p_now(i,1) = p_now(i,1)+trans10-trans01
	p_now(i,2) = p_now(i,2)-trans10+trans01
end do
end subroutine



subroutine check()
integer i, j
! Check probability is normalized
if (abs(sum(p_now) - 1.) .gt. eps) then
	write (*,*) 'Normalization Broken', sum(p_now) - 1.
	stop
end if
! Check probability is positive
do i = 1,n
do j = 1,2
	if ((p_now(i,j) .lt. 0.) .and. (abs(p_now(i,j)) .gt. eps)) &
		stop 'Negative Probability'
end do
end do
end subroutine


function drift0(x)
	real drift0, x
	drift0 = -k0*x
end function


function drift1(x)
	real drift1, x
	drift1 = -k1*(x-g)
end function


subroutine dump()
integer i

open (unit=10, file=outfile)
do i = 1, n
	write (10,*) positions(i), p_now(i,1), p_now(i,2), p_eq(i,1), p_eq(i,2)
end do
write (10,*) ''
write (10,*) ''
end subroutine

end program
