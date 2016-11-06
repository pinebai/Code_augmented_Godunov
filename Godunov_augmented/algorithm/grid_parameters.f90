module grid_and_parameters

implicit none

integer, parameter :: ighost=2, NRITER=100
real(8), parameter :: cfl=0.24d0, TOLL=0.0d0, QMAX=2.0d0
real(8), parameter :: error_data=-1.0d8

real(8) :: xl, xr
real(8):: dx, dt, rr, final_time, current_time, tol
integer :: current_step, final_step
integer :: cells_number 

real(8), dimension(:), allocatable :: x

character*12 :: boundary_type


contains


subroutine error_message

print*, 'Something is wrong here!!!'
pause

end subroutine error_message


end module grid_and_parameters