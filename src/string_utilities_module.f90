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

subroutine edit_ctrl (var, num)
! Rewrite control_file.txt
! uses: read_cmd_line, trim_zeros
character(len=25) :: buffer, label, ctrl_file, dummy
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
close(1);	close(2)
ios = rename('new.txt', ctrl_file)
end subroutine



subroutine make_filename()
end subroutine



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

end module
