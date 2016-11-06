module time_step
! Compute time step satisfying CFL condition.

use solution
! 'solution.f90'

use physics
! 'physics.f90'

implicit none


contains


!******************************************************************************************
!get dt using CFL condition number
subroutine cfldt

implicit none

integer :: i
real*8 :: rho, v, p, temp1, temp2, temp, rlamdamax

rlamdamax=0.0d0

do i=1, cells_number
 rho=u(i,1)   
 v=f_v(u(i,1),u(i,2),u(i,3))
 p=f_p(u(i,1),u(i,2),u(i,3))
 temp=dabs(v)+dsqrt(gamma*p/rho)
 if(temp>rlamdamax) rlamdamax=temp
end do

dt=dx/rlamdamax*cfl

end subroutine cfldt 


end module time_step