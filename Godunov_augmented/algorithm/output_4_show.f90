module output_show

use solution
! 'solution.f90'

implicit none

public  output_for_show


contains


subroutine output_for_show

implicit none

integer :: i

open(1,file='d:\Godunov_augmented\show\grid.dat')
write(1,'(i5)') cells_number
write(1,'(f12.6)') xl
write(1,'(f12.6)') xr
close(1)

open(2,file='d:\Godunov_augmented\show\solution.dat')
do i=1, cells_number
 write(2,'(5f15.6)'), 0.5d0*(x(i-1)+x(i)), u(i,1), v_aver(i), p_aver(i), uu(i)
end do
close(2) 

end subroutine output_for_show


end module output_show