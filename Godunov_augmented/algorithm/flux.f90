module flux

!use consant

use solution
! 'solution.f90'

use physics
! 'physics.f90'

use RiemannSolver
! 'RiemannSolver.f90'


contains
!***************************************************************************************  

! conpute f(u(i+1/2)),and f(u)=u                                                 
subroutine fluxu1_2

implicit none

integer :: i,k,IW
real*8 :: RS,US,PS
real*8 :: nu,temp,fuu1,fuu2,fuu3,fuu4
real*8,dimension(-ighost:cells_number+ighost) :: R1,U1,P1,R2,U2,P2, vstar,fuu_star
real*8,dimension(1:3) :: u_temp

real*8 :: u_star, p_star

do i=-ighost+1, cells_number+ighost 

 R1(i)=u(i,1)-d_rho(i)
 U1(i)=v_aver(i)-d_v(i)
 P1(i)=p_aver(i)-d_p(i)
      
 R2(i)=u(i,1)+d_rho(i)
 U2(i)=v_aver(i)+d_v(i)
 P2(i)=p_aver(i)+d_p(i)
!    endif
enddo

do i=0, cells_number 

 CALL RIEMANN0(R2(i),U2(i),P2(i),R1(i+1),U1(i+1),P1(i+1),u_star,p_star,IW)
 CALL RIEMANN1(0.0D0, RS,US,PS)


 vstar(i)=US

 if(RS<tol) then
  RS=tol
! 	   print*,i,"density is negative"
 endif
 if(PS<tol) then
  PS=tol
! 	   print*,i,"pressure is negative"
 endif


 flux_u(i,1)=f1(RS,RS*US,PS/(GAMMA-1.0d0)+0.5D0*RS*US**2.0D0)
 flux_u(i,2)=f2(RS,RS*US,PS/(GAMMA-1.0d0)+0.5D0*RS*US**2.0D0)
 flux_u(i,3)=f3(RS,RS*US,PS/(GAMMA-1.0d0)+0.5D0*RS*US**2.0D0)

 fuu_star(i)=f_uu(RS,US,PS)
 flux_uu(i)=f_uu(RS,US,PS)
  
end do 
   
end  subroutine fluxu1_2
!***************************************************************************************


!*************************************************************************************** 
end module flux