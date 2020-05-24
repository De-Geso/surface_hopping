! range_finder.f90
! author: Brad Friesen
! updated: 2019-11-09
! purpose: Find an appropriate lambda range and xrange for given k's and G0
! Compile: gfortran -O3 -fdefault-real-8 range_finder.f90
! Run: run_range_finder.sh

module range_finder
use parameters
use string_utilities
use functions
implicit none


contains

subroutine find_lrange ()
	lambda_max = -log(sqrt(k0/k1)*(tol/(1.0-tol))) / beta
	lambda_min = -log(sqrt(k0/k1)*(1.0/tol-1.0)) / beta
	
	write (*,*) "Lambda min:", lambda_min
	write (*,*) "Lambda max:", lambda_max
	call edit_ctrl ('lambda_min', lambda_min)
	call edit_ctrl ('lambda_max', lambda_max)
end subroutine find_lrange

subroutine find_xrange ()
	real a, b, c, fa, fb, fc, lambda
	integer i
	
	! Initial guesses for x range
	a = 0.0; b = 1000.0; c = (a + b)/2.0
	! Set lambda0 so we're in the right well
	lambda = 100.0
	
	fa = cum_prob_eq(a,lambda) - (1.0 - tol)
	fb = cum_prob_eq(b,lambda) - (1.0 - tol)
	fc = cum_prob_eq(c,lambda) - (1.0 - tol)

	! Bail out if we don't bound zero
	! write (*,*) fa, fb	
	if (fa*fb > 0.0) stop "NoBd1"
	
	! Shrink interval with golden-section search
	do i = 1, 64
		if (fa*fc < 0.0) then; b = c; fb = fc; end if
		if (fc*fb < 0.0) then; a = c; fa = fc; end if
		c = (a+b)/2.0; fc = cum_prob_eq(c, lambda) - (1.0 - tol)
		! write (*,*) b-a, c, fc
	end do
	x_max = c
	
	! Initial guesses for x range
	a = -1000.0; b = 0.0; c = (a + b)/2.0
	! Set lambda0 so we're in the left well
	lambda = -100.0
	
	fa = cum_prob_eq(a,lambda) - tol
	fb = cum_prob_eq(b,lambda) - tol
	fc = cum_prob_eq(c,lambda) - tol
	
	! Bail out if we don't bound zero
	! write (*,*) fa, fb
	if (fa*fb > 0.0) stop "NoBd2"
	
	! Shrink interval with golden-section search
	do i = 1, 64
		if (fa*fc < 0.0) then; b = c; fb = fc; end if
		if (fc*fb < 0.0) then; a = c; fa = fc; end if
		c = (a+b)/2.0; fc = cum_prob_eq(c, lambda) - tol
		! write (*,*) b-a, c, fc
	end do
	x_min = c
	
	write (*,*) 'x min:', x_min
	write (*,*) 'x max:', x_max
	call edit_ctrl ('x_min', x_min)
	call edit_ctrl ('x_max', x_max)
end subroutine find_xrange

end module range_finder
