module periodic_boundary
! Treat periodic boundary condition.

use solution
! 'solution.f90'

implicit none

public periodic_boundary_cell, periodic_boundary_HS


contains


!******************************************************************************************
! Periodic boundary half steps.
subroutine periodic_boundary_HS

implicit none

integer :: k
   
d_rho(0)=d_rho(cells_number)
d_v(0)=d_v(cells_number)
d_p(0)=d_p(cells_number) 
   
d_rho(cells_number+1)=d_rho(1)  
d_v(cells_number+1)=d_v(1)
d_p(cells_number+1)=d_p(1)  

end subroutine periodic_boundary_HS
!******************************************************************************************


!******************************************************************************************
! Periodic boundary cell-averages.
subroutine periodic_boundary_cell

implicit none

integer :: k

do k=1,3
 u(0,k)=u(cells_number,k)   
enddo 
uu(0)=uu(cells_number)  
v_aver(0)=v_aver(cells_number)
p_aver(0)=p_aver(cells_number)

do k=1,3
 u(cells_number+1,k)=u(1,k) 
enddo 
uu(cells_number+1)=uu(1)
v_aver(cells_number+1)=v_aver(1) 
p_aver(cells_number+1)=p_aver(1)    

end subroutine periodic_boundary_cell
!******************************************************************************************


end module periodic_boundary