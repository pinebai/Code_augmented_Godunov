program Blast_waves_initialize

use physics
! 'physics.f90'

use solution
! 'solution.f90'

use out_put
! 'out_put.f90'

implicit none

! Physical values
real(8) :: rhol, vl, pl, rhom, vm, pm, rhor, vr, pr

integer :: i, k, middle_cell1, middle_cell2
real*8 :: v, p, e

xl=0.0d0; xr=1.0d0

rhol=1.0d0; vl=0.0d0; pl=1000.0d0
rhom=1.0d0; vm=0.0d0; pm=1.0d-2
rhor=1.0d0; vr=0.0d0; pr=100.0d0
!middle_cell1=cells_number/10; middle_cell2=cells_number*9/10; tol=1.0d-10

print*, 'INPUT CELLS NUMBER.'
read(*,'(i5)') cells_number

middle_cell1=cells_number/10; middle_cell2=cells_number*9/10; tol=1.0d-10

allocate(x(-ighost:cells_number+ighost))

allocate(u(-ighost:cells_number+ighost,3))

allocate(uu(-ighost:cells_number+ighost)); allocate(v_aver(-ighost:cells_number+ighost)); allocate(p_aver(-ighost:cells_number+ighost))
   
! Fix grid.
dx=(xr-xl)/dfloat(cells_number)
do i=0, cells_number 
 x(i)=xl+dfloat(i)*dx
enddo

do i=1, cells_number 
 if(i.le.middle_cell1) then
  u(i,1)=rhol
  u(i,2)=f_m(rhol,vl,pl)
  u(i,3)=f_e(rhol,vl,pl)
  uu(i)=entropy_per_volume(rhol,vl,pl)
  v_aver(i)=vl
  p_aver(i)=pl 
 else
  if(i.le.middle_cell2) then
   u(i,1)=rhom
   u(i,2)=f_m(rhom,vm,pm)
   u(i,3)=f_e(rhom,vm,pm)
   uu(i)=entropy_per_volume(rhom,vm,pm)
   v_aver(i)=vm
   p_aver(i)=pm
  else
   u(i,1)=rhor
   u(i,2)=f_m(rhor,vr,pr)
   u(i,3)=f_e(rhor,vr,pr)
   uu(i)=entropy_per_volume(rhor,vr,pr)
   v_aver(i)=vr
   p_aver(i)=pr 
  end if       
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

boundary_type='reflective  '

call output_solution

!******************************************************************************************
end program Blast_waves_initialize