module physics

use grid_and_parameters
! 'grid_parameters.f90'

implicit none

real(8), parameter :: gamma=1.4d0

public  kinetic_energy, f_v, f_p, f_m, f_e, f_uu, f1, f2, f3, sgn, rminmod, entropy_per_mass, entropy_per_volume, shock_entropy_increase, gamma 


contains


!******************************************************************************************
function kinetic_energy(rho,m,e)

implicit none

real*8 :: rho, m, e, kinetic_energy

kinetic_energy=0.5d0*m**2.0d0/rho
return

end function kinetic_energy
!******************************************************************************************


!******************************************************************************************
function f_v(rho,m,e)

implicit none

real*8 :: rho, m, e, f_v

f_v=m/rho
return

end function f_v
!******************************************************************************************


!******************************************************************************************
function f_p(rho,m,e)

implicit none

real*8 :: rho, m, e, f_p

f_p=(gamma-1.0d0)*(e-kinetic_energy(rho,m,e))

return

end function f_p
!******************************************************************************************


!******************************************************************************************
function f_m(rho,v,p)

implicit none

real*8 :: rho, v, p, f_m

f_m=rho*v
return

end function  f_m
!******************************************************************************************


!******************************************************************************************
function f_e(rho,v,p)

implicit none

real*8 :: rho, v, p, f_e

f_e=p/(gamma-1.0d0)+0.5d0*rho*v*v 
return

end function f_e
!******************************************************************************************


!******************************************************************************************
function f1(rho,m,e)

implicit none

real*8 :: rho, m, e, f1

f1=m
return

end function f1
!******************************************************************************************


!******************************************************************************************
function f2(rho, m, e)

implicit none

real*8 :: rho,m,e,f2

f2=rho*f_v(rho,m,e)**2.0d0+f_p(rho,m,e)
return

end function f2
!******************************************************************************************


!******************************************************************************************
function f3(rho,m,e)

implicit none

real*8 :: rho, m, e, f3

f3=(e+f_p(rho,m,e))*f_v(rho,m,e)
return

end function f3
!******************************************************************************************


!******************************************************************************************
function sgn(x)

implicit none

real*8 :: x, sgn

if(x>0.0d0) then
   sgn=1.0d0
elseif(x<0.0d0) then
   sgn=-1.0d0
else
   sgn=0.0d0
endif

end function sgn
!******************************************************************************************


!******************************************************************************************
function rminmod(x,y)

implicit none

real*8 :: x, y, rminmod

if(x>0.0d0.and.y>0.0d0) then
 rminmod=dmin1(x,y)
else
 if(x<0.0d0.and.y<0.0d0) then
  rminmod=dmax1(x,y)
 else
  rminmod=0.0d0
 end if
endif

end function rminmod
!******************************************************************************************


!******************************************************************************************
function entropy_per_mass(rho,v,p)

implicit none

real*8 :: rho, v, p, entropy_per_mass

entropy_per_mass=dlog(p/rho**gamma)

end function entropy_per_mass
!******************************************************************************************


!******************************************************************************************
function entropy_per_volume(rho,v,p)

implicit none

real*8 :: rho, v, p, entropy_per_volume

entropy_per_volume=rho*entropy_per_mass(rho,v,p)

end function entropy_per_volume
!******************************************************************************************


!******************************************************************************************
function f_uu(rho,v,p)

implicit none

real*8 :: rho, v, p, f_uu

f_uu=v*entropy_per_volume(rho,v,p)

end function f_uu
!******************************************************************************************


!******************************************************************************************
function shock_entropy_increase(rho_l,v_l,p_l,rho_r,v_r,p_r,s)

implicit none

real*8 :: rho_l, v_l, p_l, rho_r, v_r, p_r, s, shock_entropy_increase

shock_entropy_increase=s*(entropy_per_volume(rho_r,v_r,p_r)-entropy_per_volume(rho_l,v_l,p_l))-(f_uu(rho_r,v_r,p_r)-f_uu(rho_l,v_l,p_l))

end function shock_entropy_increase
!******************************************************************************************


end module physics