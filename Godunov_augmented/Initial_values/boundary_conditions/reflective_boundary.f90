module reflective_boundary
! Treat periodic boundary condition.

use solution
! 'solution.f90'

implicit none

public reflective_boundary_cell, reflective_boundary_HS


contains


!******************************************************************************************
! Periodic boundary half steps.
subroutine reflective_boundary_HS

implicit none

integer :: k
   
d_rho(0)=d_rho(1)
d_v(0)=-d_v(1)
d_p(0)=d_p(1) 
   
d_rho(cells_number+1)=d_rho(cells_number)  
d_v(cells_number+1)=-d_v(cells_number)
d_p(cells_number+1)=d_p(cells_number)  

end subroutine reflective_boundary_HS
!******************************************************************************************


!******************************************************************************************
! Periodic boundary cell-averages.
subroutine reflective_boundary_cell

implicit none

u(0,1)=u(1,1); u(0,2)=-u(1,2); u(0,3)=u(1,3)   
uu(0)=uu(1)  
v_aver(0)=-v_aver(1)
p_aver(0)=p_aver(1)

u(cells_number+1,1)=u(cells_number,1); u(cells_number+1,2)=-u(cells_number,2); u(cells_number+1,3)=u(cells_number,3) 
uu(cells_number+1)=uu(cells_number)
v_aver(cells_number+1)=-v_aver(cells_number) 
p_aver(cells_number+1)=p_aver(cells_number)    

end subroutine reflective_boundary_cell
!******************************************************************************************


end module reflective_boundary