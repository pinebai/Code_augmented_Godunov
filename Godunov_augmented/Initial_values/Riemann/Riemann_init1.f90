program Riemann_initialize

use physics
! 'physics.f90'

use solution
! 'solution.f90'

use out_put
! 'out_put.f90'

implicit none

! Geometrical values
!real(8) :: xl, xr

! Physical values
real(8) :: rhol, vl, pl, rhor, vr, pr

integer :: i, k, middle_cell, case_number
character*3 :: if_middle_cell
real*8 :: v, p, e

xl=0.0d0; xr=1.0d0

print*, 'INPUT CASE NUMBER.'
read(*,'(i5)') case_number

print*, 'INPUT CELLS NUMBER.'
read(*,'(i5)') cells_number

!print*, 'MIDDLE CELL or NOT?'
!read(*,'(3a)') if_middle_cell

if_middle_cell='no '

allocate(x(-ighost:cells_number+ighost))

allocate(u(-ighost:cells_number+ighost,3))

allocate(uu(-ighost:cells_number+ighost)); allocate(v_aver(-ighost:cells_number+ighost)); allocate(p_aver(-ighost:cells_number+ighost))

select case(case_number)

 case(1)
  rhol=1.0d0; vl=0.0d0; pl=1.0d0
  rhor=0.125d0; vr=0.0d0; pr=0.1d0
  middle_cell=cells_number/2; tol=1.0d-12
   
 case(2) 
  rhol=1.0d0; vl=0.0d0; pl=1.0d0
  rhor=0.4263194230d0; vr=0.9274526117d0; pr=0.3031301790d0
  middle_cell=cells_number/2; tol=1.0d-12
   
 case(3)
  rhol=0.4263194230d0; vl=-0.9274526117d0; pl=0.3031301790d0
  rhor=1.0d0; vr=0.0d0; pr=1.0d0
  middle_cell=cells_number/2; tol=1.0d-12

 case(4) 
  rhol=0.2655737113d0; vl=0.9274526117d0; pl=0.3031301790d0
  rhor=0.125d0; vr=0.0d0; pr=0.1d0
  middle_cell=cells_number/2; tol=1.0d-12
   
 case(5) 
  rhol=0.4263194230d0; vl=0.9274526117d0; pl=0.3031301790d0
  rhor=0.2655737113d0; vr=0.9274526117d0; pr=0.3031301790d0
  middle_cell=cells_number/2; tol=1.0d-12
   
 case(6)
  rhol=1.0d0; vl=0.0d0; pl=1.0d0
  rhor=0.338388117766695d0; vr=1.15267997980118; pr=0.219371975711225
  middle_cell=cells_number/2; tol=1.0d-12

 case(7)
  rhor=1.0d0; vr=0.0d0; pr=1.0d0
  rhol=0.125d0; vl=0.0d0; pl=0.1d0
  middle_cell=cells_number/2; tol=1.0d-12

 case(8) 
  rhor=0.5d0; vr=0.0d0; pr=0.571d0
  rhol=0.445d0; vl=0.698d0; pl=3.528d0
  middle_cell=cells_number/2; tol=1.0d-12
   
 case(9)
  rhol=1.0d0; vl=0.75d0; pl=1.0d0
  rhor=0.125d0; vr=0.0d0; pr=0.1d0   
  middle_cell=cells_number/3; tol=1.0d-12

 case(10)
  rhol=1000.0d0; vl=0.0d0; pl=1000.0d0
  rhor=1.0d0; vr=0.0d0; pr=1.0d0
  middle_cell=cells_number/3; tol=1.0d-10

 case(11)
  rhol=7.0d0; vl=-1.0d0; pl=0.2d0
  rhor=7.0d0; vr=1.0d0; pr=0.2d0
  middle_cell=cells_number/2; tol=1.0d-12
    
 case default
    
end select

! Fix grid.
dx=(xr-xl)/dfloat(cells_number)
do i=0, cells_number 
 x(i)=xl+dfloat(i)*dx
enddo

do i=0, cells_number 
 if(i.lt.middle_cell) then
  u(i,1)=rhol
  u(i,2)=f_m(rhol,vl,pl)
  u(i,3)=f_e(rhol,vl,pl)
  uu(i)=entropy_per_volume(rhol,vl,pl)
  v_aver(i)=vl
  p_aver(i)=pl 
 else
  if(i.gt.middle_cell) then
   u(i,1)=rhor
   u(i,2)=f_m(rhor,vr,pr)
   u(i,3)=f_e(rhor,vr,pr)
   uu(i)=entropy_per_volume(rhor,vr,pr)
   v_aver(i)=vr
   p_aver(i)=pr
  else
   select case(if_middle_cell)
    case('yes')
     u(i,1)=0.5d0*(rhol+rhor)
     u(i,2)=0.5d0*(f_m(rhol,vl,pl)+f_m(rhor,vr,pr))
     u(i,3)=0.5d0*(f_e(rhol,vl,pl)+f_e(rhor,vr,pr))
     uu(i)=0.5d0*(entropy_per_volume(rhol,vl,pl)+entropy_per_volume(rhor,vr,pr))
     v_aver(i)=0.5d0*(vl+vr)
     p_aver(i)=0.5d0*(pl+pr)
    case('no ')
     u(i,1)=rhol
     u(i,2)=f_m(rhol,vl,pl)
     u(i,3)=f_e(rhol,vl,pl)
     uu(i)=entropy_per_volume(rhol,vl,pl)
     v_aver(i)=vl
     p_aver(i)=pl 
    case default; call error_message
   end select
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

boundary_type='infinitive  '

call output_solution

!******************************************************************************************
end program Riemann_initialize