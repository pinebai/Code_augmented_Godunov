program Riemann_initialize

use physics
! 'physics.f90'

use solution
! 'solution.f90'

use out_put
! 'out_put.f90'

implicit none

! Physical values
real(8) :: rhol, vl, pl, cl, edge_velocity_l, rhor, vr, pr, cr, edge_velocity_r

integer :: i, ii, case_number
real(8) :: density, momentum, energy, velocity, pressure, entropyy
real(8) :: xx, ddx, time, shift

xl=0.0d0; xr=1.0d0; tol=1.0d-12

print*, 'INPUT CASE NUMBER'
read(*,'(i5)') case_number

print*, 'INPUT CELLS NUMBER.'
read(*,'(i5)') cells_number

print*, 'INPUT TIME'
read(*,'(f20.10)') time

print*, 'INPUT POSITION SHIFT'
read(*,'(f20.10)') shift

allocate(x(-ighost:cells_number+ighost))

allocate(u(-ighost:cells_number+ighost,3))

allocate(uu(-ighost:cells_number+ighost)); allocate(v_aver(-ighost:cells_number+ighost)); allocate(p_aver(-ighost:cells_number+ighost))

select case(case_number)

 case(1)

! Left side:  
  rhol=1.0d0; vl=0.0d0; pl=2.0d0; cl=dsqrt(gamma*pl/rhol); edge_velocity_l=vl-cl
  
! Right side:
  edge_velocity_r=0.9d0 
  vr=((gamma-1)*vl+2.0d0*(cl+edge_velocity_r))/(gamma+1) 
  rhor=((rhol**gamma*(vr-edge_velocity_r)**2.0d0)/gamma/pl)**(1.0d0/(gamma-1))
  pr=pl/rhol**gamma*rhor**gamma
  
  tol=1.0d-12
   
!  time=0.1d0 
    
 case default
    
end select

! Fix grid.
dx=(xr-xl)/dfloat(cells_number)
do i=0, cells_number 
 x(i)=xl+dfloat(i)*dx
enddo

ddx=dx/1000.0d0
u=0.0d0; uu=0.0d0; v_aver=0.0d0; p_aver=0.0d0
do i=1, cells_number
 density=0.0d0; momentum=0.0d0; energy=0.0d0; velocity=0.0d0; pressure=0.0d0 !; entropy=0.0d0
 do ii=1,1000
  xx=dfloat(i-1)*dx+(dfloat(ii)-shift)*ddx
  velocity=left_wave_velocity(xx-shift,time)
  density=left_wave_density(velocity,xx-shift,time)
  pressure=left_wave_pressure(density,xx-shift,time)
  momentum=f_m(density,velocity,pressure)
  energy=f_e(density,velocity,pressure)
  entropyy=entropy_per_volume(density,velocity,pressure)
 
  u(i,1)=u(i,1)+density; u(i,2)=u(i,2)+momentum; u(i,3)=u(i,3)+energy      ! Middle rectangule formula.
  v_aver(i)=v_aver(i)+velocity; p_aver(i)=p_aver(i)+pressure; uu(i)=uu(i)+entropyy
 end do
 u(i,:)=u(i,:)/1000.0d0
 v_aver(i)=v_aver(i)/1000.0d0; p_aver(i)=p_aver(i)/1000.0d0; uu(i)=uu(i)/1000.0d0
end do

do i=1, cells_number
 velocity=f_v(u(i,1),u(i,2),u(i,3)) 
 pressure=f_p(u(i,1),u(i,2),u(i,3)) 
 entropyy=entropy_per_volume(u(i,1),velocity,pressure)
 if(p_aver(i).gt.pressure-toll) then
  if(p_aver(i).gt.pressure+tol)  call error_message 
  p_aver(i)=pressure-toll
 endif
 if(dabs(v_aver(i)-velocity).lt.tol) then
  if(v_aver(i).gt.velocity) then
   v_aver(i)=velocity+tol
  else
   v_aver(i)=velocity-tol
  end if
 end if
 if(uu(i).gt.entropyy-toll) then
  if(uu(i).gt.entropyy+tol) call error_message
  uu(i)=entropyy-toll
 end if 
end do

current_time=time
boundary_type='infinitive  '

call output_solution


contains 


function left_wave_velocity(location,time) result(a)

real(8), intent(in) :: location, time
real(8) :: a

if(location.le.time*edge_velocity_l) then
 a=vl
else 
 if(location.lt.time*edge_velocity_r) then
  a=((gamma-1.0d0)*vl+2.0d0*(cl+location/time))/(gamma+1.0d0)
 else
  a=vr
 end if
end if

end function left_wave_velocity


function left_wave_density(velocity,location,time) result(a)

real(8), intent(in) :: velocity, location, time
real(8) :: a

if(location.lt.time*edge_velocity_l) then
 a=rhol
else 
 if(location.lt.time*edge_velocity_r) then
  a=((rhol**gamma*(velocity-location/time)**2.0d0)/gamma/pl)**(1.0d0/(gamma-1.0d0))
 else
  a=rhor
 end if
end if

end function left_wave_density


function left_wave_pressure(density,location,time) result(a)

real(8), intent(in) :: density, location, time
real(8) :: a

if(location.lt.time*edge_velocity_l) then
 a=pl
else 
 if(location.lt.time*edge_velocity_r) then
  a=pl/rhol**gamma*density**gamma
 else
  a=pr
 end if
end if

end function left_wave_pressure


!******************************************************************************************
end program Riemann_initialize