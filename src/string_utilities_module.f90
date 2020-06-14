! string_utilities.f90
! Brad Friesen
! Anything useful that involves manipulating strings including file and
! command line stuff
! Compile with: gfortran -O3 -fdefault-real-8 -c string_utilities.f90
! Updated: 2020-05-10

module string_utilities
use parameters
implicit none

contains

function filelength(filename)
character(len=*) :: filename
integer i, nlines, filelength

nlines = 0
open(1, file=filename)
do
	read(1, *, end=10)
	nlines = nlines+1
end do
10 close(1)
filelength = nlines
end function



subroutine edit_ctrl (var, num)
! Rewrite control_file.txt
! uses: read_cmd_line, trim_zeros
character(len=100) :: buffer, label, ctrl_file, dummy
character(len=*) :: var
integer :: ios = 0, line = 0, col = 0
real num

var = trim(var)
ctrl_file = 'control_file.txt'

open (unit=1, file=ctrl_file)
open (unit=2, file='new.txt', status='replace')
! ios < 0 if end is encountered. ios > 0 if error detected. ios = 0 otherwise
do while (ios == 0)
	read (1, '(A)', iostat=ios) buffer
	if (ios == 0) then
		line = line + 1
		! Find first =, which is our marker
		col = scan(buffer, '=')
		! Split up label and data
		label = buffer(1:col-1)
		buffer = buffer(col+1:)
		! Copy everything to a new file
		! Trim removes trailing spaces, adjustl removes leading spaces.
		if (label == var) then
			write (*,*) 'New value for: ', trim(var), num
			write (dummy, *) num
			call trim_zeros(dummy)		! get rid of trailing zeros
			write (2,'(A)') trim(trim(label) // ' = ' // adjustl(dummy))
		! Put in blank lines
		else if (label == '') then
			write (2,*) ''
		! Rewrite lines that we don't want to change
		else
			call trim_zeros(buffer)
			write (2,'(A)') trim(trim(label) // ' = ' // adjustl(buffer))
		end if
	end if
end do
close(1)
close(2)
ios = rename('new.txt', ctrl_file)
end subroutine



function make_filename (path, s, lmbda) result (outs)
! Takes a given path and name prefix and concatenates values of k, 
! G0, lambda0, and tau if asked for to filename
! Calls: replace_text
! Updated: 2020-03-05
real, optional :: lmbda
character(*) :: path, s
character(len=100) :: outs, str
outs = s
! Write in value for k0
write (str, '(F100.5)') k0
call trim_zeros (str)
outs = trim (outs) // '_k0_' // trim (adjustl (str))
! Write in value for k1
write (str, '(F100.5)') k1
call trim_zeros (str)
outs = trim (outs) // '_k1_' // trim (adjustl (str))
! Write in our value for G0
if (G0 /= 0.0) then
	write (str, '(F100.5)') G0
	call trim_zeros (str)
	outs = trim (outs) // '_G_' // trim (adjustl (str))
end if
! Write in our value for lambda if asked for
if (present(lmbda)) then
	write (str, '(F100.5)') lmbda
	call trim_zeros (str)
	outs = trim (outs) // '_l_' // trim (adjustl (str))
end if

! add suffix, and tell us where it's going
outs = trim(outs)//'.dat'
outs = trim (path)//trim (outs)
write (*,*) 'file at: ', outs
end function



subroutine read_control()
character(*), parameter :: infile = 'control_file.txt'
character(20) buffer, label
integer ios, line, col

ios = 0
line = 0
col = 0

open (unit=10, file=infile)
! ios < 0 if end of file is encountered, ios > 0 if error has occured,
! otherwise ios = 0.
do while (ios .eq. 0)
	! '(A)' is a format specifier which says to read characters.
	read (10, '(A)', iostat=ios) buffer
	if (ios .eq. 0) then
		line = line + 1
		! Find and split label and data.
		col = scan(buffer, '=')
		label = buffer(1:col-1)
		buffer = buffer(col+1:)

		! Assign value to appropriate variable
		select case (label)
			case default
				write (*, *) 'Skipping invalid line at line: ', line
			case ('')
			! Coupling
			case ('gam')
				read (buffer, *) gam
				write (*,*) 'gam = ', gam
			! Force constants
			case ('k0')
				read (buffer, *) k0
				write (*,*) 'k0 = ', k0
			case ('k1')
				read (buffer, *) k1
				write (*,*) 'k1 = ', k1
			! Test values of x and lambda
			case ('x0')
				read (buffer, *) x0
				write (*,*) 'x0 = ', x0
			case ('lmbda0')
				read (buffer, *) lmbda0
				write (*,*) 'lmbda0 = ', lmbda0
			case ('tau')
				read (buffer, *) tau
				write (*,*) 'tau = ', tau
			case ('x_min')
				read (buffer, *) x_min
				write (*,*) 'x_min = ', x_min
			case ('x_max')
				read (buffer, *) x_max
				write (*,*) 'x_max = ', x_max
			case ('lmbda_min')
				read (buffer, *) lmbda_min
				write (*,*) 'lmbda_min = ', lmbda_min
			case ('lmbda_max')
				read (buffer, *) lmbda_max
				write (*,*) 'lmbda_max = ', lmbda_max
			! Intrinsic energy offset
			case ('G0')
				read (buffer, *) G0
				write (*,*) 'G0 = ', G0
			! Thermodynamic beta
			case ('beta')
				read (buffer, *) beta
				write (*,*) 'beta = ', beta
		end select
	end if
end do
close (10)
end subroutine



subroutine trim_zeros (string)
character (len=*) string
character (len=len(string)+2) str
integer i

str = string		! working copy of the string
if (index(str, '.') == 0) then		! if there's no decimal, add one
	str = trim(str) // '.'
end if
if (index(adjustl(str), '.') == 1) then		! if there's no number before the decimal, add a leading zero
	str = '0' // adjustl (str)
end if

do i = len_trim(str), 1, -1
	select case (str(i:i))
		case ('0')		! found a trailing 0, keep going
			cycle
		case ('.')		! found a decimal point, add a zero and exit
			str = str(1:i)
			str = trim(str) // '0'
			exit
		case default		! found a nonzero character so trim and exit
			str = str(1:i)
			str = trim(str)
			exit
	end select
end do
string = str
end subroutine trim_zeros



! Get the nth command line argument and make it nice for writing into
! control file.
function read_cmd_line (n) result (input)
	integer n
	character(len=100) :: input


	call get_command_argument(n, input)
	! bail out if we didn't supply an nth argument
	if (len_trim(input) == 0) stop "input"
	write (*,'(A, I2.0)') ' Input Number =', n
	write (*,*) 'Input Read = ', adjustl(trim(input))
end function

end module
