! edit_control.f90
! author: Brad Friesen
! updated: 2019-11-09
! purpose: Change a variable in the control file
! use: ./edit_control.out "VariableName" VariableValue

program edit_control
use parameters
use string_utilities
implicit none

character(:), allocatable :: var, str
real num

call read_control ()
! get what we want to change and what we want to change it to
var = read_cmd_line (1)
str = read_cmd_line (2)

read (str, *) num
call edit_ctrl (var, num)

end program edit_control
