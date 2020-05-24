! test_implicit.f90
! Brad Friesen
! playing with implicit PDEs

program test_implicit
use linear_algebra
implicit none

! beta = 1 is fully implicit, beta = 0.5 is Crank Nicolson
real, parameter :: beta = 1.
real, parameter :: L = 2.
real, parameter :: tau = 5.
real, parameter :: D = 10.
integer, parameter :: n = 200
real, parameter :: dt = 0.001
real, parameter :: dx = L/(n-1.)
real, parameter :: alpha = D*dt/dx**2
real, parameter :: pi = 4.*atan(1.)

real, dimension(n) :: u, v, x, exact, error
real t, norm
integer i, j



t = 0.

do i = 1,n
	x(i) = (i-1)*dx
	u(i) = f(x(i))
	v(i) = f(x(i))
end do
!u(1) = 1.0
!v(1) = 1.0
norm = sum(v)

do while (t .LE. tau)
!	exact = erfc(x/(2.*sqrt(t)))
!	error = abs(u-exact)/exact
	write (*,*) 'Normalization:', sum(v)
	call dump()
	call update()
	t = t + dt
end do

if (beta .eq. 1.) then
	call rename('fort.100', 'implicit.dat')
else if (beta .eq. 0.5) then
	call rename('fort.100', 'crank.dat')
end if



contains



subroutine update()
real, dimension(n) :: d, b
real, dimension(n-1) :: dl, du
integer i, info
real, dimension(n) :: aa, bb, cc, rr, uu
real, dimension(n,3) :: ai

dl = beta*alpha
d = -2.*beta*alpha - 1.
du = beta*alpha

do i = 2,n-1
	b(i) = -(1.-beta)*alpha*(u(i-1)-2.*u(i)+u(i+1)) - u(i)
end do

! Dirichlet boundary conditions on left
b(1) = u(1)
du(1) = 0.
d(1) = 1.

! von Neumann boundary conditions
! left
b(1) = -2.*(1.-beta)*alpha*(u(2)-u(1)) - u(1)
du(1) = 2.*beta*alpha
d(1) = -2.*beta*alpha - 1.
! right
b(n) = -2.*(1.-beta)*alpha*(u(n-1)-u(n)) - u(n)
dl(n-1) = 2.*beta*alpha
d(n) = -2.*beta*alpha - 1.

! cast things for tridag to compare
aa = 0.; aa(2:n) = dl
bb = d
cc = 0.; cc(1:n-1) = du
rr = b
ai(:,1) = aa; ai(:,2) = bb; ai(:,3) = cc;

call DGTSV (n, 1, dl, d, du, b, n, info)
if (info .NE. 0) stop info
u = b

call tridag (aa, bb, cc, rr, uu, n)
call banmul (ai, n , 1, 1, n, 3, uu, rr)
v = uu
end subroutine



subroutine test()
real, dimension(3,3) :: a
real, dimension(3) :: u, r

a(:,1) = (/0., -2., -2./)
a(:,2) = (/1., 1., 1./)
a(:,3) = (/-2., -2., 0./)
u = (/1., 1., 1./)
call banmul (a, 3, 1, 1, 3, 3, u, r)
write (*,*) r
end subroutine



subroutine improve (a, n, np, mp, b, x)
integer n, np, mp
real a(np,mp), b(n), x(n), r(n)
integer i, j
double precision sdp(n)

call banmul (a, n, 1, 1, n, 3, x, r)
r = r - b

call tridag(a(:,1), a(:,2), a(:,3), r, x, n)
x = x-r
end subroutine



function f(x)
real f, x
f = -0.5*x+1
end function



subroutine dump()
integer i

do i = 1,n
	write (100,*) x(i), u(i), v(i)
end do
write (100,*) ''
write (100,*) ''
end subroutine

end program
