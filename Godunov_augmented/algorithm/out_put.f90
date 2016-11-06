module out_put

use solution
! 'solution.f90'

implicit none


contains 


subroutine output_solution

implicit none

integer :: i

open(5,file='d:\Godunov_augmented\output\parameters.dat')
write(5,'(i5)') cells_number
write(5,'(f25.18)') tol
write(5,'(a12)') boundary_type
write(5,'(2f25.18)') xl, xr
close(5)

open(1,file='d:\Godunov_augmented\output\grid.dat')
do i=0, cells_number
 write(1,'(f25.16)') x(i)
end do
close(1)

open(2,file='d:\Godunov_augmented\output\solutions.dat')
do i=1, cells_number
 write(2,'(3f25.16)') u(i,1), u(i,2), u(i,3)
end do
close(2)

open(3,file='d:\Godunov_augmented\output\augmented.dat')
do i=1, cells_number
 write(3,'(3f25.16)') uu(i), v_aver(i), p_aver(i)
end do  
close(3)

open(4,file='d:\Godunov_augmented\output\time_step.dat')
write(4,'(f25.16)') current_time
write(4,'(i5)') current_step
close(4)

end subroutine output_solution


end module out_put