program Shu_initialize

use physics
! 'physics.f90'

use solution
! 'solution.f90'

use out_put
! 'out_put.f90'

implicit none

! Physical values
real(8) :: rhol, vl, pl
integer :: i, k, middle_cell, ii
real*8 :: v, p, e
real(8) :: xx, ddx

print*, 'INPUT CELLS NUMBER.'
read(*,'(i5)') cells_number

xl=-5.0d0; xr=5.0d0

rhol=3.857143d0; vl=2.629369d0; pl=10.33333d0

allocate(x(-ighost:cells_number+ighost))

allocate(u(-ighost:cells_number+ighost,3))

allocate(uu(-ighost:cells_number+ighost)); allocate(v_aver(-ighost:cells_number+ighost)); allocate(p_aver(-ighost:cells_number+ighost))

middle_cell=cells_number/10; tol=1.0d-10
   
! Fix grid.
dx=(xr-xl)/dfloat(cells_number); ddx=dx/10000.0d0
do i=0, cells_number 
 x(i)=xl+dfloat(i)*dx
enddo

do i=1, cells_number 
 if(i.le.middle_cell) then
  u(i,1)=rhol
  u(i,2)=f_m(rhol,vl,pl)
  u(i,3)=f_e(rhol,vl,pl)
  uu(i)=entropy_per_volume(rhol,vl,pl)
  v_aver(i)=vl
  p_aver(i)=pl 
 else
  u(i,1)=1.0d0-(dcos(5.0d0*x(i))-dcos(5.0d0*x(i-1)))/25.0d0/dx
  u(i,2)=0.0d0
  u(i,3)=2.5d0  
  do ii=0, 10000
   xx=x(i-1)+dfloat(ii)*ddx
   if(ii.ne.0.and.ii.ne.10000) then
    uu(i)=uu(i)+entropy_per_volume(1.0d0+0.2d0*dsin(5.0d0*xx),0.0d0,1.0d0)
   else
    uu(i)=uu(i)+0.5d0*entropy_per_volume(1.0d0+0.2d0*dsin(5.0d0*xx),0.0d0,1.0d0)
   end if
  end do
  uu(i)=uu(i)/10000.0d0   	   
  v_aver(i)=0.0d0
  p_aver(i)=1.0d0
      
 endif
end do

do i=1, cells_number
 v=f_v(u(i,1),u(i,2),u(i,3)) 
 p=f_p(u(i,1),u(i,2),u(i,3)) 
 e=entropy_per_volume(u(i,1),v,p)
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
 if(uu(i).gt.e-toll) then
  if(uu(i).gt.e+tol) call error_message
  uu(i)=e-toll
 end if 
end do

boundary_type='infinitive  '

call output_solution

!******************************************************************************************
end program Shu_initialize