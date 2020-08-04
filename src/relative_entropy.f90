! relative_entropy.f90 

program relative_entropy
use parameters
use string_utilities
implicit none


real eq_initial(2), eq_final(2)
real simulated(2), total, error
real x, a, b
character(:), allocatable :: outfile, infile, eq_file
integer i, ios


call init()

open (unit=10, file=eq_file)
open (unit=20, file=infile)
do while (ios .eq. 0)
	! read data
	read (10, *, iostat=ios) x, eq_initial(1), eq_initial(2), eq_final(1), eq_final(2)
	read (20, *, iostat=ios) x, simulated(1), simulated(2), a, b
	! calculate Kullback-Leibler divergence
	if (simulated(1) .gt. 0.0) then
		total = total + simulated(1)*log(simulated(1)/eq_final(1))
		error = error + log(simulated(1))+1
	end if
	if (simulated(2) .gt. 0.0) then
		total = total + simulated(2)*log(simulated(2)/eq_final(2))
		error = error + log(simulated(2))+1
	end if
end do

open (unit=30, file=outfile, position="append")
write (30, *) tau, total, abs(error)
close (30)



contains



subroutine init()
call read_control()
ios = 0
total = 0.0

eq_file = make_filename('data/distribution/', 'dist_eq')

select case (flag)
case default
stop 'Illegal flag when choosing path'
case (1)
	outfile = make_filename('data/relative_entropy/', 'rel_naive', gam=gam)
	infile = make_filename('data/distribution/', 'dist_fpe_naive', tau=tau, gam=gam)
case (2)
	outfile = make_filename('data/relative_entropy/', 'rel_fvar', gam=gam)
	infile = make_filename('data/distribution/', 'dist_fpe_fvar', tau=tau, gam=gam)
end select	
end subroutine

end program
