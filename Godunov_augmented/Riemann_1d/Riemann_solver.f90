module Riemann_solver

use fluids

implicit none

contains

subroutine riemann(ul,ur,ul_star,ur_star,wave_pattern,vacuum)

implicit none

type(phys), intent(in) :: ul, ur
type(phys), intent(out) :: ul_star, ur_star
character*11, dimension(3), intent(out) :: wave_pattern
logical, intent(out) :: vacuum

real*8 :: ml, mr
real*8 :: p_star, v_star, rhol_star, rhor_star 
real*8 :: pl, pr, rhol, rhor, vl, vr
real*8 :: sl, sr
real*8 :: al, ar
real*8 :: v_star_l, v_star_r
integer :: q, q_stop
real*8 :: ea, eb, c
real*8 :: epsilon, xxx

call tran5(rhol,vl,pl,ul) 
call tran5(rhor,vr,pr,ur)

vacuum=.false.

! vacuum check.
v_star_l=vl+integer_chara(rhol,pl,ul)
v_star_r=vr-integer_chara(rhor,pr,ur)
if(v_star_l.le.v_star_r) then
 print*, 'There is vacuum between two sides.'
 vacuum=.true.

! call phys_clon(ul_star,ul)
! call phys_clon(ur_star,ur)

! rhol_star=1.0d-10; p_star=1.0d-10
! call tran3(rhol_star,v_star_l,p_star,ul_star)
! rhor_star=1.0d-10; p_star=1.0d-10
! call tran3(rhor_star,v_star_r,p_star,ur_star)
 return 
end if

epsilon=1.0d-10
ea=1.0d-15
eb=1.0d+15

q=0
q_stop=200


do while(f(ea)*f(eb).le.0.0d0.and.dabs(ea-eb).ge.epsilon.and.q.lt.q_stop)
 c=0.5d0*(ea+eb)
 if(f(ea)*f(c).le.0.0d0) then
  eb=c
 else
  ea=c
 end if
 q=q+1
end do

p_star=c
!p_star=-800.0d0
xxx=f(p_star)

ml=phi(rhol,pl,ul,p_star)
mr=phi(rhor,pr,ur,p_star)
v_star=(pl-pr+ml*vl+mr*vr)/(ml+mr)
 
if(p_star.gt.pl)then
! The left center wave is a shock.
 sl=(rhol*vl-ml)/rhol
 rhol_star=ml/(v_star-sl)
 wave_pattern(1)='shock      '
 write(*,*)'The left center wave is a shock'
 else
!The left center wave is a rarefaction wave.
 rhol_star=density_rare(rhol,pl,ul,p_star)
 wave_pattern(1)='rarefaction'
 write(*,*)'The left center wave is a rarefaction wave'
end if

if(p_star.gt.pr)then
!The right center wave is a shock.
 sr=(rhor*vr+mr)/rhor
 rhor_star=-mr/(v_star-sr) 
 wave_pattern(3)='shock      '
 write(*,*)'The right center wave is a shock'
else
!The right center wave is a rarefaction wave.
 rhor_star=density_rare(rhor,pr,ur,p_star)
 write(*,*)'The right center wave is a rarefaction wave'
end if

wave_pattern(2)='contact   '

call phys_clon(ul_star,ul)
call phys_clon(ur_star,ur)

call tran3(rhol_star,v_star,p_star,ul_star)
call tran3(rhor_star,v_star,p_star,ur_star)


contains

 
function f(x)

implicit none

real*8 :: f,x

f=(vl-vr+pl/phi(rhol,pl,ul,x)+pr/phi(rhor,pr,ur,x))/ &
  (1.0d0/phi(rhol,pl,ul,x)+1.0d0/phi(rhor,pr,ur,x))-x

end function f

end subroutine riemann


end module Riemann_solver