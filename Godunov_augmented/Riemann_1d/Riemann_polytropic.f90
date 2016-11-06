module Riemann_polytropic_Phi
! Phi function and computation density through a rarefaction.

use polytropic

implicit none

public  phi_polytropic, density_rare_polytropic


contains


function phi_polytropic(rho,p,u,p_star)

implicit none

real(8), intent(in) :: rho, p, p_star
type(polytro), intent(in) :: u
real(8) :: phi_polytropic
real(8) :: gamma, epslon

epslon=1.0d-6
gamma=u%thermo%gamma

if(p_star/p.ge.1.0d0-epslon)then
 phi_polytropic=dsqrt(p*rho)*dsqrt((gamma+1.d0)*(p_star/p)/2.0d0 &
                +(gamma-1.d0)/2.0d0)
else
 phi_polytropic=(p_star-p)*(gamma-1.0d0)/(2.0d0*((p_star/p)**((gamma &
                -1.d0)/(2.0d0*gamma))-1.0d0)*dsqrt(gamma*p/rho))               
end if

end function phi_polytropic


function integer_chara_polytropic(rho,p,u)

implicit none

real(8), intent(in) :: rho, p
type(polytro), intent(in) :: u
real(8) :: integer_chara_polytropic
real(8) :: gamma

gamma=u%thermo%gamma

integer_chara_polytropic=2.0d0*dsqrt(gamma*p/rho)/(gamma-1.0d0)

end function integer_chara_polytropic


function density_rare_polytropic(rho_side,p_side,u_side,p_star) result(c)

implicit none

real(8), intent(in) :: p_star, rho_side, p_side
type(polytro) :: u_side
real(8) :: c
real(8) :: aa, gamma

gamma=u_side%thermo%gamma

aa=p_side/(rho_side**gamma)
c=(p_star/aa)**(1.d0/gamma)

end function density_rare_polytropic


end module Riemann_polytropic_Phi