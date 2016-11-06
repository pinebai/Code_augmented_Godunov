program main 
! Scheme for Euler system with velocity, pressure and entropy augumented.

use flux
! 'flux.f90'

use solution_reconstruction
! 'step_reconstruction.f90'

! use half_step_reconstruction

use average
! 'average.f90'

use in_put
! 'in_put.f90'

use out_put
! 'out_put.f90'

use time_step
! 'Time_step.f90'

use boundary_conditions
! 'boundary_conditions.f90'

use output_show
! 'output_4_show.f90'

implicit none

integer :: i,k,step

real*8 :: temp, rho, v, p, error_max, error_L1, temp1, temp2
real*8 :: rho_l, rho_r
real*8 :: v_l, v_r
real*8 :: p_l, p_r

real(8), allocatable, dimension(:,:) :: u0
real(8), allocatable, dimension(:) :: uu0, vv0, pp0 
! Numerical solution for temporary use.

! real(8), dimension(1:cells_number) :: xxx, yyy, zzz, www
! real(8) :: xxx, yyy

call input_solution

print*, 'INPUT FINAL STEP' 
read(*,'(i5)') final_step

print*, 'INPUT FINAL TIME'
read(*,*) final_time

dx=(xr-xl)/dfloat(cells_number)
step=0

allocate(d_rho(-ighost:cells_number+ighost)); allocate(d_v(-ighost:cells_number+ighost)); allocate(d_p(-ighost:cells_number+ighost)) 
allocate(entropy_increase(-ighost:cells_number+ighost))

allocate(d_rho_tvd(-ighost:cells_number+ighost)); allocate(d_v_tvd(-ighost:cells_number+ighost)); allocate(d_p_tvd(-ighost:cells_number+ighost)) 

allocate(flux_u(-ighost+2:cells_number+ighost-2,3)); allocate(flux_uu(-ighost+2:cells_number+ighost-2))

allocate(u0(-ighost:cells_number+ighost,3))
allocate(uu0(-ighost:cells_number+ighost)); allocate(vv0(-ighost:cells_number+ighost)); allocate(pp0(-ighost:cells_number+ighost))

do while(current_time.lt.final_time.and.step.lt.final_step)
  
 call boundary_cell_average 
  
 uu0=uu
 u0=u
 vv0=v_aver
 pp0=p_aver
  
! rho_theta0=rho_theta
  
 call cfldt
   
 if((current_time+dt).gt.final_time.and.current_time.lt.final_time) dt=final_time-current_time

 rr=dt/dx ! Set the mesh ratio.
   
 call step_reconstruction

 call boundary_HS

 call fluxu1_2 
    
 do i=1, cells_number

  rho_l=u(i,1)-d_rho(i)
  v_l=v_aver(i)-d_v(i)
  p_l=p_aver(i)-d_p(i)
   
  rho_r=u(i,1)+d_rho(i)
  v_r=v_aver(i)+d_v(i)
  p_r=p_aver(i)+d_p(i)
  
 enddo
   
 call average_vandp
  
 do i=1, cells_number
  do k=1,3
   temp=flux_u(i,k)-flux_u(i-1,k)
   u(i,k)=u(i,k)-rr*temp
  enddo
  
  temp=flux_uu(i)-flux_uu(i-1) 
  uu(i)=uu(i)-rr*temp-rr*entropy_increase(i)  
 enddo 

 call boundary_cell_average
  
 current_time=current_time+dt
 step=step+1; current_step=current_step+1
 print*, current_time, step
   
!  call output_for_show
!  pause
   
enddo 

call output_solution

call output_for_show

!***************************************************************************************
end program main