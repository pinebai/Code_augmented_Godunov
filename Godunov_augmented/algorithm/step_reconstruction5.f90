module solution_reconstruction
! Compute the half steps for the reconstruction.

use solution
! 'solution.f90'

use physics
! 'physics.f90'

implicit none

real(8), dimension(:), allocatable :: uu_memory, v_aver_memory, p_aver_memory

public  step_reconstruction
private compute_tvd_limiters, fix_limiter_reversity, compute_half_steps, conservation_trend_check, reconstruction_output,ABS_mass_HS, ABS_velocity_HS, Vaver_via_Paver_drho, &
        Paver_via_Vaver_dv, Vaver_via_drho_dv, Paver_via_drho_dv, entropy_reconstruction, v_aver_memory, p_aver_memory, uu_memory


contains


subroutine step_reconstruction

implicit none

allocate(uu_memory(1:cells_number)); allocate(v_aver_memory(1:cells_number)); allocate(p_aver_memory(1:cells_number))
v_aver_memory=v_aver(1:cells_number); p_aver_memory=p_aver(1:cells_number); uu_memory=uu(1:cells_number)

call compute_tvd_limiters
call fix_limiter_reversity
call compute_half_steps

call conservation_trend_check

call reconstruction_output

deallocate(uu_memory,v_aver_memory,p_aver_memory)

end subroutine step_reconstruction


subroutine compute_tvd_limiters
! Compute the TVD density, velocity and pressure TVD limiters.

implicit none

integer :: i

real*8 :: rho, m, e
real*8 :: rho_minus, rho_plus
real*8 :: v_minus, v, v_plus
real*8 :: p_minus, p, p_plus
   
do i=1, cells_number 
  
 rho=u(i,1)
 m=u(i,2)
 e=u(i,3) 
  
 rho_minus=u(i-1,1)
 rho_plus=u(i+1,1)
  
 v_minus=f_v(u(i-1,1),u(i-1,2),u(i-1,3)) 
 v=f_v(u(i,1),u(i,2),u(i,3)) 
 v_plus=f_v(u(i+1,1),u(i+1,2),u(i+1,3)) 
  
 p_minus=f_p(u(i-1,1),u(i-1,2),u(i-1,3)) 
 p=f_p(u(i,1),u(i,2),u(i,3)) 
 p_plus=f_p(u(i+1,1),u(i+1,2),u(i+1,3)) 

 if(p_aver(i).gt.p-toll) then
  if(p_aver(i).gt.p+tol)  call error_message 
  p_aver(i)=p-toll
 endif
 if(dabs(v_aver(i)-v).lt.tol) then
  if(v_aver(i).gt.v) then
   v_aver(i)=v+tol
  else
   v_aver(i)=v-tol
  end if
 end if
  
 d_rho_tvd(i)=rminmod(rho-rho_minus,rho_plus-rho)
 if(dabs(d_rho_tvd(i)).lt.tol) then
  if(rho_plus-rho_minus.lt.0.0d0) then
   d_rho_tvd(i)=-tol
  else
   d_rho_tvd(i)=tol
  end if
 end if

 d_v_tvd(i)=rminmod(rminmod(v_aver(i)-v_minus,v_plus-v_aver(i)),rminmod(v-v_minus,v_plus-v))
 if(dabs(d_v_tvd(i)).lt.tol) then
  if(v_plus-v_minus.lt.0.0d0) then
   d_v_tvd(i)=-tol
  else
   d_v_tvd(i)=tol
  end if
 end if

 d_p_tvd(i)=rminmod(rminmod(p_aver(i)-p_minus,p_plus-p_aver(i)),rminmod(p-p_minus,p_plus-p))
 if(dabs(d_p_tvd(i)).lt.tol) then
  if(sgn(p_plus-p_minus).lt.0.0d0) then
   d_p_tvd(i)=-tol
  else
   d_p_tvd(i)=tol
  end if
 end if 

end do

end subroutine compute_tvd_limiters


subroutine fix_limiter_reversity
! Fix the reversity of mass and velocity HS's.

implicit none

integer :: i

real*8 :: rho, m, e, v, p
logical :: yes_no

do i=1, cells_number

 rho=u(i,1)
 m=u(i,2)
 e=u(i,3) 
 v=f_v(u(i,1),u(i,2),u(i,3)) 
 p=f_p(u(i,1),u(i,2),u(i,3)) 

 if(sgn(d_rho_tvd(i)*d_v_tvd(i)*(v-v_aver(i))).lt.0.0d0) then
!  print*, '&&&&&', i
  if(v_aver(i).gt.v) then                                                        ! The TVD restriction should never be loosen.
   v_aver(i)=v-tol                                     
  else
   v_aver(i)=v+tol
  end if
!  d_v_tvd(i)=sgn(d_rho_tvd(i)*(v-v_aver(i)))*tol
   
 end if

end do

end subroutine fix_limiter_reversity


subroutine compute_half_steps
! Compute the half steps for the reconstruction.

implicit none

integer :: i
real*8 :: rho, m, e, v, p
real*8 :: drhol, drhor, drho, dpl, dpr, dv2
real(8) :: temp, temp1, temp1_p, temp2_p, temp_rho
real*8 :: d_rho_temp, d_v_temp, dp_star, drho0, dp0, paver, alpha 
real(8) :: f_mmax, f_max, f_d_rho_max, dp_max

! real(8) :: xxx

do i=1, cells_number

 rho=u(i,1)
 m=u(i,2)
 e=u(i,3) 
 v=f_v(u(i,1),u(i,2),u(i,3)) 
 p=f_p(u(i,1),u(i,2),u(i,3)) 
  
 call ABS_mass_HS(rho,v,p,v_aver(i),p_aver(i),d_rho_temp)                        ! Step1. Compute the mass HS with (7.1).

 if(dabs(d_rho_tvd(i))<d_rho_temp) then
  d_rho(i)=d_rho_tvd(i)                                                          ! Case 1.2. The mass HS computed with (7.1) does not satisfy the mass TVD restriction; therefore, 
!  flag_tvd_rho(i)=1                                                              ! the mass HS takes the TVD restricted HS and the velocisyt-average is modified accordinly. 
  call Vaver_via_Paver_drho(rho,d_rho(i),v,p,p_aver(i),v_aver(i))

!   call Paver_via_Vaver_drho(rho,d_rho(i),v,v_aver(i),p,xxx)

!****************************************************       
! 	   print*,	'Q和密度台阶修改', i,x(i)
!****************************************************
 else
  d_rho(i)=sgn(d_rho_tvd(i))*d_rho_temp                                          ! Case 1.1. The mass HS computed with (7.1) satisfies the mass TVD restriction. 
 endif
  
 call ABS_velocity_HS(rho,v,p,v_aver(i),p_aver(i),d_v_temp)                      ! Step 2. Compute the velocity HS with (7.2).

 if(dabs(d_v_tvd(i))<dabs(d_v_temp)) then
                                                                                 ! Case 2.2. The velocity HS computed with (7.2) does not satisfy the velocity TVD restriction.
  call Paver_via_Vaver_dv(rho,v,v_aver(i),d_v_tvd(i),p,temp1_p)                  ! Compute the pressure-average that make the velocity HS just satisfy the TVD restriction.
  call Paver_via_Vaver_drho(rho,d_rho_tvd(i),v,v_aver(i),p,temp2_p)              ! Compute the pressure-average that make the mass HS just satisfy the TVD restriction.
                                                                                 ! In both the computation velocity-average is not changed.  
																				 ! Modifying (increasing) the pressure-average can decrease the velocity-HS but increase the mass-HS.
  if(temp1_p<temp2_p) then

   p_aver(i)=temp1_p                                                             ! Case 2.2.1. Modify the pressure-average to make the velocity-HS just satisfy the TVD restriction and 
   d_v(i)=d_v_tvd(i)                                                             ! the mass-HS and the mass-HS is recomputed with the modified pressure-average, in which the TVD  
   call ABS_mass_HS(rho,v,p,v_aver(i),p_aver(i),d_rho_temp)                      ! restriction on mass-HS is still maintained.
   d_rho(i)=sgn(d_rho_tvd(i))*d_rho_temp
    
!****************************************************       
!   print*,	'P和速度台阶修改', i,x(i)
!****************************************************
  else
  
   d_rho(i)=d_rho_tvd(i); d_v(i)=d_v_tvd(i)                                      ! Both mass and velocity HS's reach TVD restriction.
   call Paver_via_drho_dv(rho,d_v(i),d_rho(i),p,p_aver(i)) 
   call Vaver_via_drho_dv(rho,d_v(i),d_rho(i),v,v_aver(i))
   
!****************************************************       
!    print*,'P，Q和速度台阶修改', i,x(i)
!****************************************************
  endif 
!  flag_tvd_v(i)=1
 else
  d_v(i)=sgn(d_v_tvd(i))*d_v_temp                                                ! Case 2.1. The velocity HS computed with (7.2) satisfies the velocity TVD restriction.
 endif

! Step 3. Compute pressure-HS.
     
 f_mmax=entropy_per_volume(rho,v,p)                                              ! Godunov entropy, which should always be greater than the numerical emtropy in the cell.                                                
 if(f_mmax.lt.uu(i)-tol) then
  print*, i; call error_message 
 end if
              
 f_max=entropy_per_volume(rho,v,p_aver(i))                                       ! The maxium value of the bottom-up bowl of the entropy function.
     
 dp_max=p_aver(i)*d_rho(i)/rho                                                   ! The maxium point of the entropy function with fixed mass-HS.
 f_d_rho_max=entropy_reconstruction(rho,p_aver(i),d_rho(i),dp_max)               ! The maxium value of the entropy function with fixed mass-HS.
     
 ! Case 3.1 
 if(f_max>uu(i)) then
 ! Case 3.1.1
  if(f_d_rho_max>uu(i)) then
 ! Case 3.1.1.1  
   temp=entropy_reconstruction(rho,p_aver(i),d_rho(i),d_p_tvd(i))
   if(temp<uu(i)) then
    dpl=dp_max
    dpr=d_p_tvd(i)
    call iterative_method_dp(rho,p_aver(i),uu(i),d_rho(i),dpl,dpr,dp_star) 
     
    if(((dp_star*d_p_tvd(i)).gt.0.0d0).and.(dabs(dp_star).le.dabs(d_p_tvd(i)))) then
               
     d_p(i)=dp_star	                       ! Case 3.1.1.1.1.	
	
    else
     if((dp_star*d_p_tvd(i)).le.0.0d0) then  
        
      drho0=d_rho(i)                       ! Case 3.1.1.1.2.
      dp0=0.0d0
      if(dabs(drho0)<tol) then     
       alpha=tol
      else
       alpha=dp0/drho0
      endif
       
     else
        
      drho0=d_rho(i)                       ! Case 3.1.1.1.3
      dp0=d_p_tvd(i)
      if(dabs(drho0)<tol) then     
       alpha=tol
      else
       alpha=dp0/drho0
      endif
     endif
       
     drhol=0.0d0
     drhor=drho0
     call iterative_method_drho_cs(rho,drhol,drhor,p_aver(i),alpha,uu(i),drho) 
     d_rho(i)=drho
     d_p(i)=alpha*drho
       
     call Vaver_via_Paver_drho(rho,d_rho(i),v,p,p_aver(i),v_aver(i))
     d_v(i)=rho*(v-v_aver(i))/d_rho(i)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   print*,	' Q修改 3.1.1.1.2-3', i,x(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
           
   else

    temp_rho=d_rho(i)
    temp1=entropy_reconstruction(rho,p_aver(i),d_rho_tvd(i),d_p_tvd(i))  
! 3.1.1.2.1 
    if(temp1<uu(i)) then
     drho=0.0d0                                                                 ! Case 3.1.1.2.1
     drhol=d_rho(i)
     drhor=d_rho_tvd(i)
     call iterative_method_drho(rho,drhol,drhor,p_aver(i),uu(i),d_p_tvd(i),drho)   
     d_rho(i)=drho
     d_p(i)=d_p_tvd(i)
! 3.1.1.2.2
    else
     d_rho(i)=d_rho_tvd(i)                                                      ! Case 3.1.1.2.2        
     d_p(i)=d_p_tvd(i)
     uu(i)=temp1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  print*,	' Entropy和压力台阶修改3.1.1.2.2', i,x(i)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif 
      
    call Paver_via_Vaver_drho(rho,d_rho(i),v,v_aver(i),p,paver)
    p_aver(i)=dmax1(p_aver(i),paver)
    call ABS_velocity_HS(rho,v,p,v_aver(i),p_aver(i),dv2) 
    d_v(i)=sgn(d_v_tvd(i))*dv2
     
    if(p_aver(i)>p-toll) then
     if(p_aver(i).gt.p+tol) call error_message
     p_aver(i)=p-toll
    endif
           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  print*,	' P和Q修改 3.1.1.2.1-2', i,x(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!需要处理 
   endif
! 3.1.2

  else
                       
! 3.1.2.1
   if(((dp_max*d_p_tvd(i))>0.0d0).and.(dabs(dp_max)<dabs(d_p_tvd(i)))) then      ! Case 3.1.2

    drho0=d_rho(i)                                                               ! Case 3.1.2.1
    dp0=dp_max
    if(dabs(drho0)<tol) then     
     alpha=tol
    else
     alpha=dp0/drho0
    endif		
! 3.1.2.2和3.1.2.3			
   else
    if((dp_max*d_p_tvd(i))<0.0d0) then
     drho0=d_rho(i)                                                              ! Case 3.1.2.2
     dp0=0.0d0
     if(dabs(drho0)<tol) then     
      alpha=tol
     else
      alpha=dp0/drho0
     endif
    else
     drho0=d_rho(i)                                                              ! Case 3.1.2.3  
     dp0=d_p_tvd(i)
     if(dabs(drho0)<tol) then     
      alpha=tol
     else
      alpha=dp0/drho0
     endif
    endif
   endif
    
   drhol=0.0d0
   drhor=drho0
   call iterative_method_drho_cs(rho,drhol,drhor,p_aver(i),alpha,uu(i),drho) 
   d_rho(i)=drho
   d_p(i)=alpha*drho

   call Vaver_via_Paver_drho(rho,d_rho(i),v,p,p_aver(i),v_aver(i))
   call ABS_velocity_HS(rho,v,p,v_aver(i),p_aver(i),d_v_temp)  
   d_v(i)=sgn(d_v_tvd(i))*d_v_temp
  		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! print*,' Q修改3.1.2.1-3', i,x(i)      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		      
  endif
! 3.2
 else
   
  d_rho(i)=0.0d0                                                                 ! Case 3.2.
  d_p(i)=0.0d0
   
  p_aver(i)=rho**gamma*dexp(uu(i)/rho)
!  v_aver(i)=v
    
  if(p_aver(i)>p-toll) then
   if(p_aver(i).gt.p+tol) call error_message
   p_aver(i)=p-toll
  endif
   
  call Vaver_via_Paver_drho(rho,d_rho(i),v, p,p_aver(i),v_aver(i))
  call ABS_velocity_HS(rho,v,p,v_aver(i),p_aver(i),d_v_temp)  
  d_v(i)=sgn(d_v_tvd(i))*d_v_temp

 endif

enddo

end subroutine compute_half_steps
!******************************************************************************************


!******************************************************************************************
function entropy_reconstruction(rho,p_aver,drho,dp)

implicit none

real*8, intent(in) :: rho, p_aver, drho, dp
real(8) :: entropy_reconstruction

entropy_reconstruction=0.5d0*(rho-drho)*dlog((p_aver-dp)/(rho-drho)**gamma)+0.5d0*(rho+drho)*dlog((p_aver+dp)/(rho+drho)**gamma)

end function entropy_reconstruction
!******************************************************************************************


!******************************************************************************************
function entropy_reconst_difference(rho,p_aver,uue,drho,dp) result(c)

implicit none

real*8, intent(in) :: rho, p_aver, uue, drho, dp
real(8) :: c

c=entropy_reconstruction(rho,p_aver,drho,dp)-uue

end function entropy_reconst_difference
!******************************************************************************************


!***************************************************************************************
subroutine ABS_mass_HS(rho,v,p,vaver,paver,d_rho_temp) 
! Compute the absolute value of mass HS from momentum and energy conservation.

implicit none

real*8, intent(in) :: rho, v, p, vaver, paver
real*8, intent(out) :: d_rho_temp

real*8 :: temp_p, temp_v

temp_p=p-paver
temp_v=(vaver-v)**2.0

d_rho_temp=0.5d0*rho**3.0d0*temp_v
d_rho_temp=d_rho_temp/(temp_p/(gamma-1.0d0)+0.5d0*rho*temp_v)

d_rho_temp=dsqrt(d_rho_temp)

end subroutine ABS_mass_HS


!***************************************************************************************
subroutine ABS_velocity_HS(rho,v,p,vaver,paver,d_v_temp) 
! Compute the absolute value of velocity HS from momentum and energy conservation.

implicit none

real*8, intent(in) :: rho, v, p, vaver, paver
real*8, intent(out) :: d_v_temp

d_v_temp=2.0d0*(p-paver)/((gamma-1.0d0)*rho)+(vaver-v)**2.0d0
d_v_temp=dsqrt(d_v_temp)


end subroutine ABS_velocity_HS
!***************************************************************************************


subroutine Vaver_via_Paver_drho(rho,drho,v,p,paver,vaver) 
! Modify velocity average with fixed pressure average and mass HS. The formular comes from (7.1).

implicit none

real*8, intent(in) :: rho, drho, v, p, paver
real*8, intent(out) :: vaver

real*8 :: difference

difference=2.0d0*drho**2.0d0/((gamma-1.0d0)*rho*(rho**2.0d0-drho**2.0d0))*(p-paver)

if(vaver.gt.v) then
 vaver=v+dsqrt(difference)
else
 vaver=v-dsqrt(difference)
endif

end subroutine Vaver_via_Paver_drho


subroutine Paver_via_Vaver_dv(rho,v,vaver,dv,p,paver)
! Modify pressure average with fixed velocity average and fixed velocity HS. The formular comes from (7.2).

implicit none

real(8), intent(in) :: rho, v, vaver, dv, p
real(8), intent(out) :: paver

paver=p-0.5d0*(gamma-1.0d0)*rho*(dv**2.0d0-(vaver-v)**2.0d0)  

end subroutine Paver_via_Vaver_dv


subroutine Paver_via_Vaver_drho(rho,drho,v,vaver,p,paver)
! Modify pressure average with fixed velocity average and fixed mass HS. The formular comes from (7.1).

implicit none

real(8), intent(in) :: rho, drho, v, vaver, p
real(8), intent(out) :: paver

paver=p-0.5d0*(gamma-1.0d0)*rho*(vaver-v)**2.0d0*(rho**2.0d0-drho**2.0d0)/(drho**2.0d0)

end subroutine Paver_via_Vaver_drho


!***************************************************************************************
subroutine Vaver_via_drho_dv(rho,dv,drho,v,vaver) 

implicit none

real*8, intent(in) :: rho, dv, drho, v
real(8), intent(inout) :: vaver

vaver=v-drho*dv/rho

end subroutine Vaver_via_drho_dv
!***************************************************************************************


!***************************************************************************************
subroutine Paver_via_drho_dv(rho,dv,drho,p,paver) 

implicit none

real*8, intent(in) :: rho, dv, drho, p
real(8), intent(inout) :: paver

paver=p-dv**2.0d0*(gamma-1.0d0)*(rho**2.0d0-drho**2.0d0)/(2.0d0*rho)

end subroutine Paver_via_drho_dv
!***************************************************************************************


!***************************************************************************************
subroutine iterative_method_drho_cs(rho,drhol,drhor,paver,alpha,uue,drho) 
! Iteration solves both mass and pressure HS's on 'CS' curve, with mass HS as the parameter.

implicit none

real(8), intent(in) :: rho, drhol, drhor, paver, alpha, uue
real(8), intent(out) :: drho

real*8 :: f_valuel, f_valuer, f_valuem
real*8 :: xl, xm, xr, xm0, error
integer :: iteration_step

xl=drhol
xr=drhor
error=1.0d3
iteration_step=0

f_valuel=entropy_reconst_difference(rho,paver,uue,xl,alpha*xl)
f_valuer=entropy_reconst_difference(rho,paver,uue,xr,alpha*xr)

xm=0.5d0*(xl+xr)

if(f_valuer<=0.0d0) then
 if(f_valuel*f_valuer<=0.0d0) then
  do while(error>tol)
   xm0=xm
   f_valuel=entropy_reconst_difference(rho,paver,uue,xl,alpha*xl)
   f_valuer=entropy_reconst_difference(rho,paver,uue,xr,alpha*xr)
   f_valuem=entropy_reconst_difference(rho,paver,uue,xm,alpha*xm)
   if((f_valuel*f_valuem)>0.0d0.and.(f_valuer*f_valuem)>0.0d0) then
    call error_message		  
   endif
   if((f_valuel*f_valuem)<0.0d0.and.(f_valuer*f_valuem)<0.0d0) then
    call error_message		  
   endif
   if((f_valuel*f_valuem)<=0.0d0) then
    xr=xm
   else
    xl=xm
   endif
   xm=0.5d0*(xl+xr)
   error=dabs(xm-xm0)
   iteration_step=iteration_step+1
	     !print*,error,iterate_step
  enddo
  drho=xm
 else
  drho=0.0d0
!       drho=drhor
 endif
else
 drho=drhor
endif

end subroutine iterative_method_drho_cs
!***************************************************************************************


!***************************************************************************************
subroutine iterative_method_drho(rho,drhol,drhor,paver,uue,dp,drho) 
! Iteration solves mass HS with fixed pressure HS.


implicit none

real*8, intent(in) :: rho, drhol, drhor, paver, uue, dp
real(8), intent(out) :: drho

real*8 :: f_valuel, f_valuem, f_valuer
real*8 :: xl, xm, xr, xm0, error
integer :: iterate_step

xl=drhol
xr=drhor
error=1.0d3
iterate_step=0

f_valuel=entropy_reconst_difference(rho,paver,uue,xl,dp)
f_valuer=entropy_reconst_difference(rho,paver,uue,xr,dp)
xm=0.5d0*(xl+xr)

if(f_valuer<=0.0d0) then
 if(f_valuel*f_valuer<=0.0d0) then
  do while(error>tol)
   xm0=xm
   f_valuel=entropy_reconst_difference(rho,paver,uue,xl,dp)
   f_valuer=entropy_reconst_difference(rho,paver,uue,xr,dp)
    
   f_valuem=entropy_reconst_difference(rho,paver,uue,xm,dp)
   if((f_valuel*f_valuem)>0.0d0.and.(f_valuer*f_valuem)>0.0d0) then
    print*,'It is wrong!'
    stop		  
   endif
   if((f_valuel*f_valuem)<0.0d0.and.(f_valuer*f_valuem)<0.0d0) then
    print*,'It is wrong!'
    stop		  
   endif
   if((f_valuel*f_valuem)<=0.0d0) then
    xr=xm
   else
    xl=xm
   endif
   xm=0.5d0*(xl+xr)
   error=dabs(xm-xm0)
   iterate_step=iterate_step+1
	     !print*,error,iterate_step
  enddo
  drho=xm
 else
!       print*,"cuowu"
  drho=0.0d0
!       drho=drhor
 endif
else
   drho=drhor
endif

end subroutine iterative_method_drho
!***************************************************************************************


!***************************************************************************************
subroutine iterative_method_dp(rho,paver,uue,drho,dpl,dpr,dp) 

implicit none

real*8, intent(in) :: rho, paver, uue, drho, dpl, dpr
real(8), intent(out) :: dp

real*8 :: f_valuel, f_valuem, f_valuer
real*8 :: xl, xm, xr, xm0, error
integer :: iterate_step

xl=dpl; xr=dpr; error=1.0d3; iterate_step=0

f_valuel=entropy_reconst_difference(rho,paver,uue,drho,xl)
f_valuer=entropy_reconst_difference(rho,paver,uue,drho,xr)

xm=0.5d0*(xl+xr)

if(f_valuel*f_valuer<=0.0d0) then
 do while(error>tol)
  xm0=xm
  f_valuel=entropy_reconst_difference(rho,paver,uue,drho,xl)
  f_valuer=entropy_reconst_difference(rho,paver,uue,drho,xr)
   
  f_valuem=entropy_reconst_difference(rho,paver,uue,drho,xm)
  if((f_valuel*f_valuem)>0.0d0.and.(f_valuer*f_valuem)>0.0d0) call error_message
  if((f_valuel*f_valuem)<0.0d0.and.(f_valuer*f_valuem)<0.0d0) call error_message

  if((f_valuel*f_valuem)<=0.0d0) then
   xr=xm
  else
   xl=xm
  endif
  xm=0.5d0*(xl+xr)
  error=dabs(xm-xm0)
  iterate_step=iterate_step+1
     !print*,error,iterate_step
 enddo
 dp=xm
else
 call error_message !dp=0.0d0
endif

! continue

end subroutine iterative_method_dp


subroutine energy_relation_check(rho,drho,v,vaver,p,paver,yes_no)
! Check the energy relation (3.1) in "scratch3RR" or (4.20) in "Euler_threestep".

real(8), intent(in) :: rho, drho, v, vaver, p, paver
logical, intent(out) :: yes_no

real(8) :: xxx

yes_no=.false.

xxx=drho**2.0d0*((p-paver)/(gamma-1.0d0)+0.5d0*rho*(v-vaver)**2.0d0)-0.5d0*rho**3.0d0*(v-vaver)**2.0d0
if(dabs(xxx).lt.toll*1.0d2) yes_no=.true.

end subroutine energy_relation_check


subroutine conservation_trend_check
! Check the conservation and entropy non-decrease after the reconstruction.

implicit none

integer :: i
real(8) :: xxx, yyy, zzz
real(8) :: v, p, ss

! Conservation and entropy increase check.
do i=1, cells_number
   
 xxx=0.5d0*(f_m(u(i,1)-d_rho(i),v_aver(i)-d_v(i),p_aver(i)-d_p(i))+f_m(u(i,1)+d_rho(i),v_aver(i)+d_v(i),p_aver(i)+d_p(i)))
 yyy=0.5d0*(f_e(u(i,1)-d_rho(i),v_aver(i)-d_v(i),p_aver(i)-d_p(i))+f_e(u(i,1)+d_rho(i),v_aver(i)+d_v(i),p_aver(i)+d_p(i)))
 zzz=0.5d0*(entropy_per_volume(u(i,1)-d_rho(i),v_aver(i)-d_v(i),p_aver(i)-d_p(i))+entropy_per_volume(u(i,1)+d_rho(i),v_aver(i)+d_v(i),p_aver(i)+d_p(i)))
 if(dabs(xxx-u(i,2)).gt.1.0d-6) then
  print*, i, '   momentum'; call error_message
 end if
 if(dabs(yyy-u(i,3)).gt.1.0d-6) then
  print*, i, '   energy'; call error_message
 end if
 if(zzz-uu(i).lt.-1.0d-6) then
  print*, i, '   entropy'; call error_message
 end if
end do

! Trend of modification of velocity, pressure and entropy cell-averages check
do i=1, cells_number
 v=f_v(u(i,1),u(i,2),u(i,3))
 if(dabs(v_aver(i)-v).gt.dabs(v_aver_memory(i)-v)+tol*1.0d1) call error_message
 p=f_p(u(i,1),u(i,2),u(i,3))
 if(p_aver(i).gt.p+tol.or.p_aver(i).lt.p_aver_memory(i)-tol*1.0d1) call error_message
 ss=entropy_per_volume(u(i,1),v,p)
 if(uu(i).gt.ss+tol.or.uu(i).lt.uu_memory(i)-tol*1.0d1) call error_message
end do

end subroutine conservation_trend_check


subroutine reconstruction_output
! Output the information of the reconstructed solution for visualization.

integer :: i

open(1,file='d:\Godunov_augmented\show\reconstruction.dat')
do i=1, cells_number
 write(1,'(7f12.6)') 0.5d0*(x(i-1)+x(i)), u(i,1),d_rho(i),v_aver(i),d_v(i),p_aver(i),d_p(i)
end do
close(1)

end subroutine reconstruction_output


!***************************************************************************************
end module solution_reconstruction