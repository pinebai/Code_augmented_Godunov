program Sinewave_initialize

use physics
! 'physics.f90'

use solution
! 'solution.f90'

use out_put
! 'out_put.f90'

implicit none


! Physical values
real(8) :: pi, density, velocity, pressure, entropyy
integer :: i, ii, case_number
real(8) :: xx, ddx

xl=0.0d0; xr=1.0d0; tol=1.0d-12

print*, 'INPUT CASE NUMBER'
read(*,'(i5)') case_number

print*, 'INPUT CELLS NUMBER.'
read(*,'(i5)') cells_number

pi=4.0d0*datan(1.0d0)

allocate(x(-ighost:cells_number+ighost))

allocate(u(-ighost:cells_number+ighost,3))

allocate(uu(-ighost:cells_number+ighost)); allocate(v_aver(-ighost:cells_number+ighost)); allocate(p_aver(-ighost:cells_number+ighost))

dx=(xr-xl)/dfloat(cells_number)
do i=0, cells_number 
 x(i)=xl+dfloat(i)*dx
enddo

ddx=dx/10000.0d0

select case(case_number)

 case(1)
    
  do i=1, cells_number
   u(i,1)=1.0d0-0.1d0*(dcos(2.0d0*pi*x(i))-dcos(2.0d0*pi*x(i-1)))/dx/pi
   v_aver(i)=1.0d0; p_aver(i)=1.0d0
   u(i,2)=1.0d0-0.1d0*(dcos(2.0d0*pi*x(i))-dcos(2.0d0*pi*x(i-1)))/dx/pi
   u(i,3)=0.5d0*(1.0d0-0.1d0*(dcos(2.0d0*pi*x(i))-dcos(2.0d0*pi*x(i-1)))/dx/pi)+1.0d0/(gamma-1.0d0)
   uu(i)=0.0d0
   do ii=0, 10000
    xx=(x(i-1)+dfloat(ii)*ddx)*2.0d0*pi
    if(ii.ne.0.and.ii.ne.10000) then
     uu(i)=uu(i)+entropy_per_volume(1.0d0+0.2d0*dsin(xx),1.0d0,1.0d0)
    else
     uu(i)=uu(i)+0.5d0*entropy_per_volume(1.0d0+0.2d0*dsin(xx),1.0d0,1.0d0)
    end if
   end do
   uu(i)=uu(i)/10000.0d0   		 	 
  end do
   
 case(2)

  do i=1, cells_number
   u(i,1)=1.0d0-0.1d0*(dcos(2.0d0*pi*x(i))-dcos(2.0d0*pi*x(i-1)))/dx/pi
   v_aver(i)=0.0d0
   p_aver(i)=1.0d0-0.1d0*(dcos(2.0d0*pi*x(i))-dcos(2.0d0*pi*x(i-1)))/dx/pi
   u(i,2)=0.0d0
   u(i,3)=p_aver(i)/(gamma-1.0d0)
   uu(i)=0.0d0
   do ii=0, 10000
    xx=(x(i-1)+dfloat(ii)*ddx)*2.0d0*pi
    if(ii.ne.0.and.ii.ne.10000) then
     uu(i)=uu(i)+entropy_per_volume(1.0d0+0.2d0*dsin(xx),0.0d0,1.0d0+0.2d0*dsin(xx))
    else
     uu(i)=uu(i)+0.5d0*entropy_per_volume(1.0d0+0.2d0*dsin(xx),0.0d0,1.0d0+0.2d0*dsin(xx))
    end if
   end do
   uu(i)=uu(i)/10000.0d0   		 	 
  end do

 case(3)

  do i=1, cells_number
   u(i,1)=1.0d0-0.1d0*(dcos(2.0d0*pi*(x(i)+0.25d0))-dcos(2.0d0*pi*(x(i-1)+0.25d0)))/dx/pi
   v_aver(i)=0.0d0
   p_aver(i)=1.0d0-0.1d0*(dcos(2.0d0*pi*(x(i)+0.25d0))-dcos(2.0d0*pi*(x(i-1)+0.25d0)))/dx/pi
   u(i,2)=0.0d0
   u(i,3)=p_aver(i)/(gamma-1.0d0)
   uu(i)=0.0d0
   do ii=0, 10000
    xx=(x(i-1)+0.25d0+dfloat(ii)*ddx)*2.0d0*pi
    if(ii.ne.0.and.ii.ne.10000) then
     uu(i)=uu(i)+entropy_per_volume(1.0d0+0.2d0*dsin(xx),0.0d0,1.0d0+0.2d0*dsin(xx))
    else
     uu(i)=uu(i)+0.5d0*entropy_per_volume(1.0d0+0.2d0*dsin(xx),0.0d0,1.0d0+0.2d0*dsin(xx))
    end if
   end do
   uu(i)=uu(i)/10000.0d0   		 	 
  end do

 case(4)
    
  do i=1, cells_number
   u(i,1)=1.0d0-0.1d0*(dcos(2.0d0*pi*(x(i)+0.5d0))-dcos(2.0d0*pi*(x(i-1)+0.5d0)))/dx/pi
   v_aver(i)=1.0d0; p_aver(i)=1.0d0
   u(i,2)=1.0d0-0.1d0*(dcos(2.0d0*pi*(x(i)+0.5d0))-dcos(2.0d0*pi*(x(i-1)+0.5d0)))/dx/pi
   u(i,3)=0.5d0*(1.0d0-0.1d0*(dcos(2.0d0*pi*(x(i)+0.5d0))-dcos(2.0d0*pi*(x(i-1)+0.5d0)))/dx/pi)+1.0d0/(gamma-1.0d0)
   uu(i)=0.0d0
   do ii=0, 10000
    xx=(x(i-1)+0.5d0+dfloat(ii)*ddx)*2.0d0*pi
    if(ii.ne.0.and.ii.ne.10000) then
     uu(i)=uu(i)+entropy_per_volume(1.0d0+0.2d0*dsin(xx),1.0d0,1.0d0)
    else
     uu(i)=uu(i)+0.5d0*entropy_per_volume(1.0d0+0.2d0*dsin(xx),1.0d0,1.0d0)
    end if
   end do
   uu(i)=uu(i)/10000.0d0   		 	 
  end do

 case default; call error_message
  
end select

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

current_time=0.0d0
boundary_type='periodic    '

call output_solution

!******************************************************************************************
end program Sinewave_initialize