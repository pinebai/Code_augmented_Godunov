module in_put

use solution
! 'solution.f90'

implicit none


contains 


subroutine input_solution

implicit none

integer :: i

open(5,file='d:\Godunov_augmented\output\parameters.dat')
read(5,'(i5)') cells_number
read(5,'(f25.18)') tol
read(5,'(a12)') boundary_type
read(5,'(2f25.18)') xl, xr
close(5)

allocate(x(-ighost:cells_number+ighost))

allocate(u(-ighost:cells_number+ighost,3))

allocate(uu(-ighost:cells_number+ighost)); allocate(v_aver(-ighost:cells_number+ighost)); allocate(p_aver(-ighost:cells_number+ighost))

open(1,file='d:\Godunov_augmented\input\grid.dat')
do i=0, cells_number
 read(1,'(f25.16)') x(i)
end do
close(1)

open(2,file='d:\Godunov_augmented\input\solutions.dat')
do i=1, cells_number
 read(2,'(3f25.16)') u(i,1), u(i,2), u(i,3)
end do  
close(2)

open(3,file='d:\Godunov_augmented\input\augmented.dat')
do i=1, cells_number
 read(3,'(3f25.16)')  uu(i), v_aver(i), p_aver(i)
end do  
close(3)

open(4,file='d:\Godunov_augmented\input\time_step.dat')
read(4,'(f25.16)') current_time
read(4,'(i5)') current_step
close(4)

end subroutine input_solution


end module in_put