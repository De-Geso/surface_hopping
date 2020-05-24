! fpe.f90
! Brad Friesen
! Solve the Fokker Planck equation
! Updated: 2020-05-13

program fpe
use functions
use parameters
use string_utilities
use linear_algebra
implicit none

! theta = 1.0 for implicit, theta = 0.5 for Crank-Nicolson
real, parameter :: theta = 1.
! transfer constant. Controls how agressively probability is transfered
! between potentials. Higher is more agressive. 0.<=gam<=1.
real, parameter :: gam = 0.1
! number of times to transfer probability between potentials before
! going to the next time step
integer, parameter :: transfers = 1

real, parameter :: xx_min = -2.0
real, parameter :: xx_max = 3.0
integer, parameter :: n = 1000
real, parameter :: dt = 0.001
real, parameter :: dx = (xx_max-xx_min)/(n-1.0)
integer, parameter :: check_step = int(1.0/dt)
real, parameter :: eps = 1.e-12
real, dimension(n) :: positions
real, dimension(n, 2) :: prob, p_now, p_last, p_last_ref, &
	pot_at_pos, drift, diff
real z
integer i, j

! TESTING
real :: lmbda(n) = 0.

call init()

! Populate arrays
do i = 1, n
	pot_at_pos(i, 1) = V0(positions(i),lmbda(1))
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
!prob(:,1) = 2./(3.*n)
!prob(:,2) = 1./(3.*n)

p_now = prob
write (*,*) sum(p_now) - 1.


call steady_state()



contains



subroutine init()
! Initialize variables
integer i, j
call read_control()
prob = 0.0
p_now = 0.0
p_last = 0.0
p_last_ref = 0.0
pot_at_pos = 0.0
diff = D
z = 0.0

do i = 1, n
	positions(i) = xx_min + dx*(i-1.)
	drift(i,1) = drift0(positions(i))
	drift(i,2) = drift1(positions(i))
end do
end subroutine


subroutine steady_state()
real tot_var_dist
integer continue_condition, step_counter
integer i, j


call dump()
tot_var_dist = 0.0
continue_condition = 1
step_counter = 0

do while (continue_condition .eq. 1)
	! Save previous distribution
	p_last = p_now
	! Advance probability one timestep
	call update_prob()
!	do i = 1, transfers
!		call update_transfer(lmbda(1))
!	end do
	
	if (step_counter .ge. check_step) then
		! output the current distribution
		call dump()
		! check if the distribution changed from last check
		tot_var_dist = 0.5*sum(abs(p_last_ref-p_now))
		write (*, *) tot_var_dist, eps
		! check normalization
		if (abs(sum(p_now) - 1.) .gt. eps) then
			write (*,*) 'Normalization Broken', sum(p_now) - 1.
!			stop
		end if
		! check probability is positive
		do i = 1,n
		do j = 1,2
			if ((p_now(i,j) .lt. 0.) .and. (abs(p_now(i,j)) .gt. eps)) stop 'Negative Probability'
		end do
		end do

		if (tot_var_dist .lt. eps) then
			continue_condition = 0
		else
			tot_var_dist = 0.0
			step_counter = 0
			p_last_ref = p_now
		end if
	end if
	step_counter = step_counter + 1
end do
end subroutine



subroutine update_prob()
real dl(n-1,2), d(n,2), du(n-1,2), b(n,2), a(n,2), c(n,2), u(n,2), r(n,2)
real, dimension(n,2) :: aa, bb, cc, rr, uu
integer info(2)
integer i, j

do i = 1, n-1
	du(i,:) = -theta*dt*(diff(i+1,:)/(dx*dx) - drift(i+1,:)/(2.*dx))
	dl(i,:) = -theta*dt*(diff(i,:)/(dx*dx) + drift(i,:)/(2.*dx))
end do
d = 1. + theta*dt*2.*diff/(dx*dx)

! b is overwritten, and we still maybe care about p_last, so feed sgtsv
! a dummy
do i = 2, n-1
	b(i,:) = p_last(i,:) + dt*(1.-theta)*(		&
		-(drift(i+1,:)*p_last(i+1,:) - drift(i-1,:)*p_last(i-1,:))/(2.*dx)		&
		+(diff(i+1,:)*p_last(i+1,:)-2.*diff(i,:)*p_last(i,:)+diff(i-1,:)*p_last(i-1,:))/(dx*dx))
end do

! dirichlet boundary conditions
b(1,:) = p_last(1,:)
du(1,:) = 0.
d(1,:) = 1.

b(n,:) = p_last(n,:)
d(n,:) = 1.
dl(n-1,:) = 0.

! von Neumann boundary conditions
b(1,:) = p_last(1,:) + (1.-theta)*dt*2.*(diff(2,:)*p_last(2,:)-diff(1,:)*p_last(1,:))/(dx*dx)
d(1,:) = 1. + theta*dt*2.*diff(1,:)/(dx*dx)
du(1,:) = -theta*dt*2.*diff(2,:)/(dx*dx)

b(n,:) = p_last(n,:) + (1.-theta)*dt*2.*(diff(n-1,:)*p_last(n-1,:)-diff(n,:)*p_last(n,:))/(dx*dx)
d(n,:) = 1. + theta*dt*2.*diff(n,:)/(dx*dx)
dl(n-1,:) = -theta*dt*2.*diff(n,:)/(dx*dx)

!! von Neumann boundary conditions 2
!! left
!b(1,:) = p_last(1,:) + (1.-theta)*dt*(		&
!	(-drift(2,:)*p_last(2,:)+drift(1,:)*p_last(1,:))/(2.*dx)		&
!	+(diff(2,:)*p_last(2,:)-2.*diff(1,:)*p_last(1,:)		&
!		+diff(1,:)*(p_last(1,:)-drift(1,:)*dx/diff(1,:)*p_last(1,:)))/(dx*dx))
!d(1,:) = 1. + theta*dt*2.*diff(1,:)
!du(1,:) = -theta*dt*2.*diff(2,:)

aa = 0.; aa(2:n,:) = dl
bb = d
cc = 0.; cc(1:n-1,:) = du
rr = b

do i = 1,2
!	call tridag(aa(:,i), bb(:,i), cc(:,i), rr(:,i), uu(:,i), n)
	call dgtsv(n, 1, dl(:,i), d(:,i), du(:,i), b(:,i), n, info(i))
	if (info(i) .NE. 0) stop info(i)
end do

p_now = b
!p_now = uu

end subroutine



subroutine update_transfer(lmbda_now)
real de, lmbda_now, p_at_pos
integer i

do i = 1, n
	de = V1(positions(i))-V0(positions(i),lmbda_now)
	p_at_pos = p_now(i,1) + p_now(i,2)
	
	p_now(i,1) = (1.-fermi(de))*p_at_pos
	p_now(i,2) = (fermi(de))*p_at_pos
end do
end subroutine


function drift0(x)
	real drift0, x
	drift0 = -k0*x
	drift0 = 0.
end function


function drift1(x)
	real drift1, x
	drift1 = -k1*(x-g)
	drift1 = 0.
end function


subroutine dump()
integer i, j

do j = 1, 2
	do i = 1, n
		write (j,*) positions(i), p_now(i,j)
	end do
	write (j,*) ''
	write (j,*) ''
end do
end subroutine

end program
